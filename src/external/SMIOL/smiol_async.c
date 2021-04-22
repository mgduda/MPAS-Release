#include <stdio.h>
#include <stdlib.h>
#include "smiol_async.h"


/********************************************************************************
 *
 * SMIOL_async_init
 *
 * Short, one-line description.
 *
 * Detailed description.
 *
 ********************************************************************************/
int SMIOL_async_init(struct SMIOL_context *context)
{
	int ierr;
	pthread_mutexattr_t mutexattr;
	pthread_condattr_t condattr;


	/*
	 * Mutex setup
	 */
	context->buf_mutex = malloc(sizeof(pthread_mutex_t));

	ierr = pthread_mutexattr_init(&mutexattr);
	if (ierr) {
		fprintf(stderr, "Error: pthread_mutexattr_init: %i\n", ierr);
		return 1;
	}

	ierr = pthread_mutex_init(context->buf_mutex, (const pthread_mutexattr_t *)&mutexattr);
	if (ierr) {
		fprintf(stderr, "Error: pthread_mutex_init: %i\n", ierr);
		return 1;
	}

	ierr = pthread_mutexattr_destroy(&mutexattr);
	if (ierr) {
		fprintf(stderr, "Error: pthread_mutexattr_destroy: %i\n", ierr);
		return 1;
	}


	/*
	 * Condition variable setup
	 */
	context->buf_cond = malloc(sizeof(pthread_cond_t));

	ierr = pthread_condattr_init(&condattr);
	if (ierr) {
		fprintf(stderr, "Error: pthread_condattr_init: %i\n", ierr);
		return 1;
	}

	ierr = pthread_cond_init(context->buf_cond, (const pthread_condattr_t *)&condattr);
	if (ierr) {
		fprintf(stderr, "Error: pthread_cond_init: %i\n", ierr);
		return 1;
	}

	ierr = pthread_condattr_destroy(&condattr);
	if (ierr) {
		fprintf(stderr, "Error: pthread_condattr_destroy: %i\n", ierr);
		return 1;
	}

        context->buf_usage = 0;
        context->max_buf_usage = ((size_t)2048 * (size_t)1024 * (size_t)1024 - 1);
        context->n_bufs = 0;

	/*
	 * Asynchronous queue initialization
	 */
	context->write_queue = malloc(sizeof(struct SMIOL_async_queue));
	*(context->write_queue) = SMIOL_ASYNC_QUEUE_INITIALIZER;

	context->write_priority_queue = malloc(sizeof(heap_t));
	*(context->write_priority_queue) = HEAP_INITIALIZER;
	heap_init(context->write_priority_queue, (size_t)1000000);

	context->pending_queue = malloc(sizeof(struct SMIOL_async_queue));
	*(context->pending_queue) = SMIOL_ASYNC_QUEUE_INITIALIZER;

	/*
	 * Ticket lock initialization
	 */
	context->write_queue_lock = malloc(sizeof(struct SMIOL_async_ticketlock));
	*(context->write_queue_lock) = SMIOL_ASYNC_TICKETLOCK_INITIALIZER;
	SMIOL_async_ticketlock_create(context->write_queue_lock);


	/*
	 * Asynchronous writer thread initialization
	 */
	context->writer = NULL;


	/*
	 * Asynchronous status initialization
	 */
	context->active = 0;


	return 0;
}


/********************************************************************************
 *
 * SMIOL_async_finalize
 *
 * Short, one-line description.
 *
 * Detailed description.
 *
 ********************************************************************************/
int SMIOL_async_finalize(struct SMIOL_context *context)
{
	int ierr;

	ierr = pthread_mutex_destroy(context->buf_mutex);
	if (ierr) {
		fprintf(stderr, "Error: pthread_mutex_destroy: %i\n", ierr);
		return 1;
	}

	free(context->buf_mutex);
	context->buf_mutex = NULL;

	ierr = pthread_cond_destroy(context->buf_cond);
	if (ierr) {
		fprintf(stderr, "Error: pthread_cond_destroy: %i\n", ierr);
		return 1;
	}

	free(context->buf_cond);
	context->buf_cond = NULL;

	/*
	 * Free queue
	 */
/* TO DO: could check here that queues are empty */
	free(context->write_queue);
	heap_free(context->write_priority_queue);
	free(context->write_priority_queue);
	free(context->pending_queue);

	/*
	 * Free ticket lock
	 */
	SMIOL_async_ticketlock_free(context->write_queue_lock);
	free(context->write_queue_lock);

	return 0;
}


/********************************************************************************
 *
 * SMIOL_async_queue_add
 *
 * Short, one-line description.
 *
 * Detailed description.
 *
 ********************************************************************************/
void SMIOL_async_queue_add(struct SMIOL_async_queue *queue, struct SMIOL_async_buffer *b)
{
	if (!queue->head && !queue->tail) {
		queue->head = b;
		queue->tail = b;
	} else if (queue->head && queue->tail) {
		queue->head->next = b;
		queue->head = b;
	} else {
		fprintf(stderr, "List error: only one of head or tail was associated!\n");
	}
}


/********************************************************************************
 *
 * SMIOL_async_queue_empty
 *
 * Short, one-line description.
 *
 * Detailed description.
 *
 ********************************************************************************/
int SMIOL_async_queue_empty(struct SMIOL_async_queue *queue)
{
	if (!queue->tail) {
		return 1;
	}

	return 0;
}


/********************************************************************************
 *
 * SMIOL_async_queue_remove
 *
 * Short, one-line description.
 *
 * Detailed description.
 *
 ********************************************************************************/
struct SMIOL_async_buffer *SMIOL_async_queue_remove(struct SMIOL_async_queue *queue)
{
	struct SMIOL_async_buffer *b;

	if (!queue->tail) {
		b = NULL;
	} else {
		b = queue->tail;
		queue->tail = queue->tail->next;
		if (!queue->tail) {
			queue->head = NULL;
		}
	}

	return b;
}


/********************************************************************************
 *
 * SMIOL_async_ticketlock_create
 *
 * https://stackoverflow.com/questions/12685112/pthreads-thread-starvation-caused-by-quick-re-locking
 *
 * Detailed description.
 *
 ********************************************************************************/
void SMIOL_async_ticketlock_create(struct SMIOL_async_ticketlock *lock)
{
	int ierr;
        pthread_mutexattr_t mutexattr;
	pthread_condattr_t condattr;

        ierr = pthread_mutexattr_init(&mutexattr);
        if (ierr) {
                fprintf(stderr, "Error: pthread_mutexattr_init: %i\n", ierr);
                return;
        }

        lock->mutex = malloc(sizeof(pthread_mutex_t));
        ierr = pthread_mutex_init(lock->mutex, (const pthread_mutexattr_t *)&mutexattr);
        if (ierr) {
                fprintf(stderr, "Error: pthread_mutex_init: %i\n", ierr);
                return;
        }

        ierr = pthread_mutexattr_destroy(&mutexattr);
        if (ierr) {
                fprintf(stderr, "Error: pthread_mutexattr_destroy: %i\n", ierr);
                return;
        }

	/*
	 * Condition variable setup
	 */
	lock->cond = malloc(sizeof(pthread_cond_t));

	ierr = pthread_condattr_init(&condattr);
	if (ierr) {
		fprintf(stderr, "Error: pthread_condattr_init: %i\n", ierr);
		return;
	}

	ierr = pthread_cond_init(lock->cond, (const pthread_condattr_t *)&condattr);
	if (ierr) {
		fprintf(stderr, "Error: pthread_cond_init: %i\n", ierr);
		return;
	}

	ierr = pthread_condattr_destroy(&condattr);
	if (ierr) {
		fprintf(stderr, "Error: pthread_condattr_destroy: %i\n", ierr);
		return;
	}

	lock->queue_head = 0;
	lock->queue_tail = 0;
}


/********************************************************************************
 *
 * SMIOL_async_ticketlock_free
 *
 * https://stackoverflow.com/questions/12685112/pthreads-thread-starvation-caused-by-quick-re-locking
 *
 * Detailed description.
 *
 ********************************************************************************/
void SMIOL_async_ticketlock_free(struct SMIOL_async_ticketlock *lock)
{
	int ierr;

	/*
	 * Free mutex
	 */
        ierr = pthread_mutex_destroy(lock->mutex);
        if (ierr) {
                fprintf(stderr, "Error: pthread_mutex_destroy: %i\n", ierr);
		return;
	}

        free(lock->mutex);

	/*
	 * Free condition variable
	 */
        ierr = pthread_cond_destroy(lock->cond);
        if (ierr) {
                fprintf(stderr, "Error: pthread_cond_destroy: %i\n", ierr);
		return;
	}

        free(lock->cond);
}


/********************************************************************************
 *
 * SMIOL_async_ticket_lock
 *
 * https://stackoverflow.com/questions/12685112/pthreads-thread-starvation-caused-by-quick-re-locking
 *
 * Detailed description.
 *
 ********************************************************************************/
void SMIOL_async_ticket_lock(struct SMIOL_async_ticketlock *lock)
{
	unsigned long queue_me;

	pthread_mutex_lock(lock->mutex);
	queue_me = lock->queue_tail++;
	while (queue_me != lock->queue_head) {
		pthread_cond_wait(lock->cond, lock->mutex);
	}
	pthread_mutex_unlock(lock->mutex);
}


/********************************************************************************
 *
 * SMIOL_async_ticket_unlock
 *
 * https://stackoverflow.com/questions/12685112/pthreads-thread-starvation-caused-by-quick-re-locking
 *
 * Detailed description.
 *
 ********************************************************************************/
void SMIOL_async_ticket_unlock(struct SMIOL_async_ticketlock *lock)
{
	pthread_mutex_lock(lock->mutex);
	lock->queue_head++;
	pthread_cond_broadcast(lock->cond);
	pthread_mutex_unlock(lock->mutex);
}


/********************************************************************************
 *
 * SMIOL_async_launch_thread
 *
 * Short, one-line description.
 *
 * Detailed description.
 *
 ********************************************************************************/
void SMIOL_async_launch_thread(pthread_t **thread, void *(*func)(void *), void *arg)
{
	int ierr;
	pthread_attr_t attr;

	*thread = malloc(sizeof(pthread_t));

	ierr = pthread_attr_init(&attr);
	if (ierr) {
		fprintf(stderr, "Error: pthread_attr_init: %i\n", ierr);
	}

	ierr = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	if (ierr) {
		fprintf(stderr, "Error: pthread_attr_setdetachedstate: %i\n", ierr);
	}

	ierr = pthread_create(*thread, (const pthread_attr_t *)&attr, func, arg);
	if (ierr) {
		fprintf(stderr, "Error: pthread_create: %i\n", ierr);
	}

	ierr = pthread_attr_destroy(&attr);
	if (ierr) {
		fprintf(stderr, "Error: pthread_attr_destroy: %i\n", ierr);
	}
}


/********************************************************************************
 *
 * SMIOL_async_join_thread
 *
 * Short, one-line description.
 *
 * Detailed description.
 *
 ********************************************************************************/
void SMIOL_async_join_thread(pthread_t **thread)
{
	int ierr;

	if (*thread != NULL) {
		ierr = pthread_join(**thread, NULL);
		if (ierr) {
			fprintf(stderr, "Error: pthread_create: %i\n", ierr);
		}

		free(*thread);
		*thread = NULL;
	}
}

