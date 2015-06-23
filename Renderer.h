#include <pebble.h>

#define PI 3.14159265
#define TwoPI 6.28318530
#define FastSQRT 0x5f3759df  

#define Vector_Up ((Vector){(float)0.0,(float)1.0,(float)0.0})
#define identity_matrix ((Matrix){(float)1.0,(float)0.0,(float)0.0,(float)0.0,(float)0.0,(float)1.0,(float)0.0,(float)0.0,(float)0.0,(float)0.0,(float)1.0,(float)0.0,(float)0.0,(float)0.0,(float)0.0,(float)1.0})

typedef struct Vector
{
	float x;
	float y;
	float z;
}Vector;

typedef struct Vector2
{
	float x;
	float y;
}Vector2;

typedef struct {
	float m11, m12, m13, m14;
	float m21, m22, m23, m24;
	float m31, m32, m33, m34;
	float m41, m42, m43, m44;
} Matrix;

typedef struct {
	Vector position;
	Vector target;
}Camera;

typedef struct {
	uint8_t a; uint8_t b; uint8_t c;
}Triangle;

typedef struct {
	Vector position;
	Vector rotation;

	//vertices
	uint8_t vertices_count;
	Vector* vertices;

	//triangles
	uint8_t triangles_count;
	Triangle* triangles;

	Vector2* projection;
}Mesh;

//
// FAST MATH ---------------------------------------------------------
//

static float fastsin(float x){
	float sin;
  
  if (x < -PI)
    x += TwoPI;
  
  if (x > PI)
    x -= TwoPI;
  
  if (x < 0) {
    sin = 1.2732395 * x + .40528473 * x * x;
    
    if (sin < 0)
      sin = .225 * (sin *-sin - sin) + sin;
    else
      sin = .225 * (sin * sin - sin) + sin;
    
  } else {
      sin = 1.2732395 * x - .40528473 * x * x;
    
      if (sin < 0)
        sin = .225 * (sin *-sin - sin) + sin;
      else
        sin = .225 * (sin * sin - sin) + sin;
  }
  
  return sin;
}

inline static float fastcos(float x){
	return fastsin(x + 1.57079632);
}

inline int abs(const int i){
	return (i < 0) ? -i : i;
}

inline int absf(const float i){
	return (i < 0.0) ? -i : i;
}

inline float fastsqrt(const float x){
	union{
		float f;
		int i;
	}u;
  
	u.f = x;
	u.i = FastSQRT - (u.i >> 1);

	return x * u.f * (1.5f - 0.5f * x * u.f * u.f);
}

//
// VECTOR ---------------------------------------------------------
//

inline float length(Vector* v){
	const float x = v->x;
	const float y = v->y;
	const float z = v->z;

	return fastsqrt(x*x+y*y+z*z);
}

inline void normalize(Vector* v){
	const float scale = 1.0 / length(v);
	v->x *= scale;
	v->y *= scale;
	v->z *= scale;
}

inline void set_vector(Vector* v, float x, float y, float z){
	v->x=x,v->y=y,v->z=z;
}

//c = a - b
inline void sub(Vector* a, Vector* b, Vector* c){
	c->x = a->x - b->x;
	c->y = a->y - b->y;
	c->z = a->z - b->z;
}

//c = a * b 
inline void cross(Vector* a, Vector* b, Vector* c){
	c->x = a->y * b->z - a->z * b->y;
	c->y = a->z * b->x - a->x * b->z;
	c->z = a->x * b->y - a->y * b->x;
}

inline float dot(Vector* a, Vector* b){
	return a->x * b->x + a->y * b->y + a->z * b->z;
}

//
// MATRIX ---------------------------------------------------------
//

inline void set_matrix(Matrix* m,
					   float m11, float m12, float m13, float m14,
					   float m21, float m22, float m23, float m24,
					   float m31, float m32, float m33, float m34,
					   float m41, float m42, float m43, float m44){

	m->m11=m11; m->m12=m12; m->m13=m13; m->m14=m14;
	m->m21=m21; m->m22=m22; m->m23=m23; m->m24=m24;
	m->m31=m31; m->m32=m32; m->m33=m33; m->m34=m34;
	m->m41=m41; m->m42=m42; m->m43=m43; m->m44=m44;
}

void lookat_LH(Vector* eye, Vector* target, Vector* up, Matrix* result){
	static Vector vX; static Vector vY;  static Vector vZ;

	sub(target,eye,&vZ);
	normalize(&vZ);

	cross(up,&vZ,&vX);
	normalize(&vX);

	cross(&vZ,&vX,&vY);
	normalize(&vY);

	float x = -dot(&vX,eye);  float y = -dot(&vY,eye);  float z = -dot(&vZ,eye);

	set_matrix(result,
			   vX.x, vY.x, vZ.x,0,
			   vX.y, vY.y, vZ.y,0,
			   vX.z, vY.z, vZ.z,0,
			   x, y, z, 1);
}

void perspective_LH(float fov, float aspect, float znear, float zfar, Matrix* m){
	const float tan = 2.56;
  
  m->m11 = tan / aspect;
  m->m12 = m->m13 = m->m14 = 0.0;
  
  m->m22 = tan;
  m->m21 = m->m23 = m->m24 = 0.0;
  
  m->m31 = m->m32 = 0.0;
  m->m33 = -zfar / (znear - zfar);
  m->m34 = 1.0;
  
  m->m41 = m->m42 = m->m43 = 0.0;
  m->m44 = (znear * zfar) / (znear- zfar);
}

void matrix_multiply(Matrix* a, Matrix* b, Matrix* c){
	c->m11 = a->m11 * b->m11 + a->m12 * b->m21 + a->m13 * b->m31 + a->m14 * b->m41;
	c->m12 = a->m11 * b->m12 + a->m12 * b->m22 + a->m13 * b->m32 + a->m14 * b->m42;
	c->m13 = a->m11 * b->m13 + a->m12 * b->m23 + a->m13 * b->m33 + a->m14 * b->m43;
	c->m14 = a->m11 * b->m14 + a->m12 * b->m24 + a->m13 * b->m34 + a->m14 * b->m44;

	c->m21 = a->m21 * b->m11 + a->m22 * b->m21 + a->m23 * b->m31 + a->m24 * b->m41;
	c->m22 = a->m21 * b->m12 + a->m22 * b->m22 + a->m23 * b->m32 + a->m24 * b->m42;
	c->m23 = a->m21 * b->m13 + a->m22 * b->m23 + a->m23 * b->m33 + a->m24 * b->m43;
	c->m24 = a->m21 * b->m14 + a->m22 * b->m24 + a->m23 * b->m34 + a->m24 * b->m44;

	c->m31 = a->m31 * b->m11 + a->m32 * b->m21 + a->m33 * b->m31 + a->m34 * b->m41;
	c->m32 = a->m31 * b->m12 + a->m32 * b->m22 + a->m33 * b->m32 + a->m34 * b->m42;
	c->m33 = a->m31 * b->m13 + a->m32 * b->m23 + a->m33 * b->m33 + a->m34 * b->m43;
	c->m34 = a->m31 * b->m14 + a->m32 * b->m24 + a->m33 * b->m34 + a->m34 * b->m44;

	c->m41 = a->m41 * b->m11 + a->m42 * b->m21 + a->m43 * b->m31 + a->m44 * b->m41;
	c->m42 = a->m41 * b->m12 + a->m42 * b->m22 + a->m43 * b->m32 + a->m44 * b->m42;
	c->m43 = a->m41 * b->m13 + a->m42 * b->m23 + a->m43 * b->m33 + a->m44 * b->m43;
	c->m44 = a->m41 * b->m14 + a->m42 * b->m24 + a->m43 * b->m34 + a->m44 * b->m44;
}

void yaw_pitch_roll(float yaw, float pitch, float roll, Matrix* m){
	const float x = (fastcos(yaw*.5) * fastsin(pitch*.5) * fastcos(roll*.5)) + (fastsin(yaw*.5) * fastcos(pitch*.5) * fastsin(roll*.5));
	const float y = (fastsin(yaw*.5) * fastcos(pitch*.5) * fastcos(roll*.5)) - (fastcos(yaw*.5) * fastsin(pitch*.5) * fastsin(roll*.5));
	const float z = (fastcos(yaw*.5) * fastcos(pitch*.5) * fastsin(roll*.5)) - (fastsin(yaw*.5) * fastsin(pitch*.5) * fastcos(roll*.5));
	const float w = (fastcos(yaw*.5) * fastcos(pitch*.5) * fastcos(roll*.5)) + (fastsin(yaw*.5) * fastsin(pitch*.5) * fastsin(roll*.5));

	m->m11 = 1.0 - (2.0 * (y*y + z*z));  
  m->m12 = 2.0 * (x*y + z*w);  
  m->m13 = 2.0 * (z*x - y*w);  
  m->m14 = 0.0;
  
	m->m21 = 2.0 * (x*y - z*w);  
  m->m22 = 1.0 - (2.0 * (z*z + x*x));  
  m->m23 = 2.0 * (y*z + x*w);  
  m->m24 = 0.0;
  
	m->m31 = 2.0 * (z*x + y*w);  
  m->m32 = 2.0 * (y*z - x*w);   
  m->m33 = 1.0 - (2.0*(y*y + x*x));  
  m->m34 = 0.0;
  
	m->m41 = 0.0;  m->m42 = 0.0;  m->m43 = 0.0;  m->m44 = 1.0;
}

void translation(float x, float y, float z, Matrix* m){
	m->m11 = 1.0;  m->m12 = 0.0;  m->m13 = 0.0;  m->m14 = 0.0;
	m->m21 = 0.0;  m->m22 = 1.0;  m->m23 = 0.0;  m->m24 = 0.0;
	m->m31 = 0.0;  m->m32 = 0.0;  m->m33 = 1.0;  m->m34 = 0.0;
	m->m41 = x;    m->m42 = y;    m->m43 = z;    m->m44 = 1.0;
}

void transform_coords(Vector* v, Matrix* m, Vector2* r){
	r->x = v->x * m->m11 + v->y * m->m21 + v->z * m->m31 + m->m41;
	r->y = v->x * m->m12 + v->y * m->m22 + v->z * m->m32 + m->m42;

	r->x /= (v->x * m->m14 + v->y * m->m24 + v->z * m->m34 + m->m44);
	r->y /= (v->x * m->m14 + v->y * m->m24 + v->z * m->m34 + m->m44);
}

//
// MESH
//

inline void set_triangle(Triangle* t, uint8_t a, uint8_t b, uint8_t c){
	t->a=a; t->b=b; t->c=c;
}

static void free_mesh(Mesh* mesh){
	free(mesh->vertices);
	free(mesh->triangles);
	free(mesh->projection);
	free(mesh);
}

//
// DRAW
//

//calculate view and projection matrices
inline void calculate_predraw(Camera* cam, int w, int h, Matrix* viewMatrix, Matrix* projectionMatrix){
	lookat_LH(&cam->position, &cam->target, &Vector_Up, viewMatrix);
	perspective_LH(.78, ((float)w)/((float)h), .01, 1.0, projectionMatrix);
}

inline void draw_line(GContext* ctx, int x0, int y0, int x1, int y1){
	//Bresenham algorythm found on wikipedia
	const int dx = abs(x1 - x0);
	const int dy = abs(y1 - y0);
	const int sx = (x0 < x1) ? 1 : -1;
	const int sy = (y0 < y1) ? 1 : -1;
	int err = dx - dy;
	int e2;

	while(1){
		graphics_draw_pixel(ctx, (GPoint) {x0,y0});

		if ((x0 == x1) && (y0 == y1)) 
			break;

		e2 = err << 1;

		if(e2 > -dy){
			err -= dy;
			x0 += sx;
		}

		if(e2 < dx){
			err += dx;
			y0 += sy;
		}
	}
}

inline void fix_mesh_rotation(Mesh* mesh){
	while(mesh->rotation.x < -PI)
		mesh->rotation.x += TwoPI;
	while(mesh->rotation.x > -PI)
		mesh->rotation.x -= TwoPI;

	while(mesh->rotation.y < -PI)
		mesh->rotation.y += TwoPI;
	while(mesh->rotation.y > -PI)
		mesh->rotation.y -= TwoPI;

	while(mesh->rotation.z < -PI)
		mesh->rotation.z += TwoPI;
	while(mesh->rotation.z > -PI)
		mesh->rotation.z -= TwoPI;
}

void draw_mesh(GContext* ctx,
		  Mesh* mesh,
		  Matrix* viewMatrix, Matrix* projectionMatrix,
		  int w, int h,
		  Matrix* worldMatrix, Matrix* translateMatrix, Matrix* transformMatrix){

	uint8_t i;
	const uint8_t vertices_count = mesh->vertices_count;
	const uint8_t triangles_count = mesh->triangles_count;

	uint8_t a, b, c;
	uint8_t ax, ay, bx, by, cx, cy;

	fix_mesh_rotation(mesh);

	Vector2 projection;
  
  yaw_pitch_roll(mesh->rotation.x, mesh->rotation.y, mesh->rotation.z, worldMatrix);
	translation(mesh->position.x, mesh->position.y, mesh->position.z, translateMatrix);
		
	matrix_multiply(worldMatrix,translateMatrix,transformMatrix);
	matrix_multiply(transformMatrix,viewMatrix,transformMatrix);
	matrix_multiply(transformMatrix,projectionMatrix,transformMatrix);
  
	for(i=0; i<vertices_count;i++){
		transform_coords(&mesh->vertices[i], transformMatrix, &projection);

		mesh->projection[i].x = (projection.x * w + w * .5);
		mesh->projection[i].y = (projection.y * h + h * .5);
	}

	for(i=0; i<triangles_count;i++){
		a = mesh->triangles[i].a;
		b = mesh->triangles[i].b;
		c = mesh->triangles[i].c;

		ax = (mesh->projection[a].x);
		ay = (mesh->projection[a].y);
		bx = (mesh->projection[b].x);
		by = (mesh->projection[b].y);
		cx = (mesh->projection[c].x);
		cy = (mesh->projection[c].y);

		draw_line(ctx, ax, ay, bx, by);
		draw_line(ctx, ax, ay, cx, cy);
		draw_line(ctx, bx, by, cx, cy);
	}
}
