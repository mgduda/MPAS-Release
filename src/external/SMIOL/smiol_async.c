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
	context->mutex = malloc(sizeof(pthread_mutex_t));

	ierr = pthread_mutexattr_init(&mutexattr);
	if (ierr) {
		fprintf(stderr, "Error: pthread_mutexattr_init: %i\n", ierr);
		return 1;
	}

	ierr = pthread_mutex_init(context->mutex, (const pthread_mutexattr_t *)&mutexattr);
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
	context->cond = malloc(sizeof(pthread_cond_t));

	ierr = pthread_condattr_init(&condattr);
	if (ierr) {
		fprintf(stderr, "Error: pthread_condattr_init: %i\n", ierr);
		return 1;
	}

	ierr = pthread_cond_init(context->cond, (const pthread_condattr_t *)&condattr);
	if (ierr) {
		fprintf(stderr, "Error: pthread_cond_init: %i\n", ierr);
		return 1;
	}

	ierr = pthread_condattr_destroy(&condattr);
	if (ierr) {
		fprintf(stderr, "Error: pthread_condattr_destroy: %i\n", ierr);
		return 1;
	}

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

	ierr = pthread_mutex_destroy(context->mutex);
	if (ierr) {
		fprintf(stderr, "Error: pthread_mutex_destroy: %i\n", ierr);
		return 1;
	}

	free(context->mutex);
	context->mutex = NULL;

	ierr = pthread_cond_destroy(context->cond);
	if (ierr) {
		fprintf(stderr, "Error: pthread_cond_destroy: %i\n", ierr);
		return 1;
	}

	free(context->cond);
	context->cond = NULL;

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

