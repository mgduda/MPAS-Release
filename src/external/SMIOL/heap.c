#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "heap.h"

struct heap_item {
	void *item;
	long key;
	size_t order;
};

#define HEAP_ITEM_INITIALIZER (struct heap_item){NULL, 0, 0}

#if 0
int main(void)
{
	size_t i;
	heap_t h = HEAP_INITIALIZER;
	struct buffer *temp;
	struct buffer *min;
	long key;
	long min_key = 0;

	srand(time(NULL));

	heap_init(&h, 100000000);

	printf("Heap size = %i\n", (int)heap_size(h));

	for (i = 0; i < 100000000; i++) {
		temp = malloc(sizeof(struct buffer));
		temp->name = malloc(sizeof(char) * 32);
		snprintf(temp->name, 32, "uReconstructZonal"); 
		key = rand() % 100000;
		if (i % 3 == 0) key = -key;
		temp->deadline = key;
//		printf("insert %f\n", (double)temp->deadline);

		if (i == 0 || key < min_key) {
			min_key = key;
		}

		heap_insert(&h, (void *)temp, key);
	}

	printf("min inserted key is %li\n", min_key);

	for (i = 0; i < 100000000; i++) {
		min = heap_get_min(&h);
//		printf("min key = %f\n", (double)min->deadline);
		if (min->deadline < min_key) {
			printf("ERROR: Keys are not increasing!\n");
		}
		min_key = min->deadline;

		free(min);
	}

	heap_free(&h);
	
	return 0;
}
#endif


void heap_init(heap_t *heap, size_t heap_size)
{
	size_t i;

	heap->max_size = heap_size;
	heap->cur_size = 0;
	heap->next_order = 0;
	heap->heap = malloc(sizeof(struct heap_item) * heap->max_size);

	for (i = 0; i < heap->max_size; i++) {
		heap->heap[i] = HEAP_ITEM_INITIALIZER;
	}
}


void heap_free(heap_t *heap)
{
	free(heap->heap);
	heap->max_size = 0;
	heap->cur_size = 0;
	heap->next_order = 0;
}


void heap_insert(heap_t *heap, void *item, long key)
{
	size_t cursor;
	size_t parent;
	struct heap_item temp;

	if (heap->cur_size == heap->max_size) {
		fprintf(stderr, "Error: heap is full!\n");
		return;
	}

	cursor = heap->cur_size;
	parent = (cursor - 1) / 2;
	heap->heap[cursor].item = item;
	heap->heap[cursor].key = key;
	heap->heap[cursor].order = heap->next_order++;
	heap->cur_size++;

	while (cursor != 0 &&
	       (heap->heap[cursor].key < heap->heap[parent].key ||
	        (heap->heap[cursor].key == heap->heap[parent].key && heap->heap[cursor].order < heap->heap[parent].order))) {

		/* Swap current node with parent */
		temp = heap->heap[cursor];
		heap->heap[cursor] = heap->heap[parent];
		heap->heap[parent] = temp;
		
		/* Parent becomes new current node */
		cursor = parent;
		parent = (cursor - 1) / 2;
	}

	
}


void *heap_get_min(heap_t *heap)
{
	size_t cursor;
	size_t left, right, min;
	struct heap_item temp;
	void *item;

	if (heap->cur_size == 0) {
//		fprintf(stderr, "Error: heap is empty!\n");
		return NULL;
	}

	item = heap->heap[0].item;

	heap->cur_size--;
	heap->heap[0] = heap->heap[heap->cur_size];

	/* Could clear memory from node that was just moved to root */

	/* If size is now zero, we could return */

	min = 0;
	do {
		cursor = min;
		left = 2 * cursor + 1;
		right = 2 * cursor + 2;
		if (left < heap->cur_size) {
			if (heap->heap[left].key < heap->heap[min].key ||
			    (heap->heap[left].key == heap->heap[min].key && heap->heap[left].order < heap->heap[min].order)) {
				min = left;
			}

			/* For a complete tree, right can only exist if left exists */
			if (right < heap->cur_size) {
				if (heap->heap[right].key < heap->heap[min].key ||
				    (heap->heap[right].key == heap->heap[min].key && heap->heap[right].order < heap->heap[min].order)) {
					min = right;
				}
			}
		}

		if (min != cursor) {
			/* Swap smaller child with parent */
			temp = heap->heap[cursor];
			heap->heap[cursor] = heap->heap[min];
			heap->heap[min] = temp;
		}
	} while (min != cursor);

	return item;
}


size_t heap_size(heap_t heap)
{
	return heap.cur_size;
}
