#include <pebble.h>
#include "Renderer.h"

#define TIMER 33

static Window *window;
static Layer *layer;

int frames = 0;

Mesh *cube;
Camera *cam;

Matrix *viewMatrix;
Matrix *projectionMatrix;
Matrix *transformMatrix;
Matrix *worldMatrix;
Matrix *translateMatrix;

// UPDATE LOOP
static void update(){
  cube->rotation.x -= .05;
  cube->rotation.y += .05;
}

// DRAW
static void draw(Layer *l, GContext* ctx){
  const GRect bounds = layer_get_bounds(l);
  const uint8_t h = bounds.size.h;
  const uint8_t w = bounds.size.w;

  calculate_predraw(cam, w, h, viewMatrix, projectionMatrix);

  draw_mesh(ctx, cube, viewMatrix, projectionMatrix, w, h, worldMatrix, translateMatrix, transformMatrix);

  ++frames;  
}

static void create_cube(Mesh* mesh, float length){
  mesh->position = ((Vector){(float)0.0,(float)0.0,(float)0.0});
  mesh->rotation = ((Vector){(float)0.0,(float)0.0,(float)0.0});

  const float n = .5 * length;
  mesh->vertices_count = 8;
  mesh->vertices = malloc(8 * sizeof(Vector));
  mesh->projection = malloc(8 * sizeof(Vector2));

  set_vector(&mesh->vertices[0], -n, n, n);
  set_vector(&mesh->vertices[1], n, n, n);
  set_vector(&mesh->vertices[2], -n, -n, n);
  set_vector(&mesh->vertices[3], n, -n, n);

  set_vector(&mesh->vertices[4], -n, n, -n);
  set_vector(&mesh->vertices[5], n, n, -n);
  set_vector(&mesh->vertices[6], n, -n, -n);
  set_vector(&mesh->vertices[7], -n, -n, -n);

  mesh->triangles_count = 12;
  mesh->triangles = malloc(12 * sizeof(Triangle));

  set_triangle(&mesh->triangles[0], 0, 1, 2);
  set_triangle(&mesh->triangles[1], 1, 2, 3);
  set_triangle(&mesh->triangles[2], 1, 3, 6);
  set_triangle(&mesh->triangles[3], 1, 5, 6);
  set_triangle(&mesh->triangles[4], 0, 1, 4);
  set_triangle(&mesh->triangles[5], 1, 4, 5);

  set_triangle(&mesh->triangles[6], 2, 3, 7);
  set_triangle(&mesh->triangles[7], 3, 6, 7);
  set_triangle(&mesh->triangles[8], 0, 2, 7);
  set_triangle(&mesh->triangles[9], 0, 4, 7);
  set_triangle(&mesh->triangles[10], 4, 5, 6);
  set_triangle(&mesh->triangles[11], 4, 6, 7);
}

// INIT

static void init_visuals(){
  viewMatrix = malloc(sizeof(Matrix));
  projectionMatrix = malloc(sizeof(Matrix));
  transformMatrix = malloc(sizeof(Matrix));
  worldMatrix = malloc(sizeof(Matrix));
  translateMatrix = malloc(sizeof(Matrix));

  cam = malloc(sizeof(Camera));
  cam->target = ((Vector){(float)0.0,(float)0.0,(float)0.0});
  set_vector(&(cam->position), 0, 0, 15);

  cube = malloc(sizeof(Mesh));
  create_cube(cube, 2);
  cube->position.y = 0;

}

static void destroy(){
  free(cam);
  free_mesh(cube);
  free(viewMatrix);
  free(projectionMatrix);
  free(transformMatrix);
  free(worldMatrix);
  free(translateMatrix);
}

static void window_load(Window *window) {
  Layer *window_layer = window_get_root_layer(window);

  GRect bounds = layer_get_bounds(window_layer);
  layer = layer_create(bounds);
  layer_set_update_proc(layer, draw);
  layer_add_child(window_layer, layer);
}

static void window_unload(Window *window) {
  layer_destroy(layer);
}

void handle_timer(void *data){
  app_timer_register(TIMER, handle_timer, NULL);
  update();
  layer_mark_dirty(layer);
}

static void handle_tick_s(struct tm *tick, TimeUnits units_changed){
  APP_LOG(APP_LOG_LEVEL_DEBUG, "FPS: %d", frames);
  frames = 0;
}

static void init(void) {
  window = window_create();
  window_set_window_handlers(window, (WindowHandlers) {
    .load = window_load,
    .unload = window_unload,
  });
  window_set_fullscreen(window, true);
  window_stack_push(window, false);

  tick_timer_service_subscribe(SECOND_UNIT, handle_tick_s);
}

static void deinit(void) {
  window_destroy(window);
}

int main(void) {
  init_visuals();
  init();

  app_timer_register(TIMER, handle_timer, NULL);
  app_event_loop();
  deinit();
  destroy();
}
