#ifndef SMIOL_HEAP_H
#define SMIOL_HEAP_H

typedef struct {
	size_t max_size;
	size_t cur_size;
	size_t next_order;
	struct heap_item *heap;
} heap_t;

#define HEAP_INITIALIZER (heap_t){0, 0, 0, NULL}


void heap_init(heap_t *heap, size_t heap_size);
void heap_free(heap_t *heap);
void heap_insert(heap_t *heap, void *item, long key);
void *heap_get_min(heap_t *heap);
size_t heap_size(heap_t heap);

#endif
