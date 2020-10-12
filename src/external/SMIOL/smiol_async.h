#ifndef SMIOL_ASYNC_H
#define SMIOL_ASYNC_H

#include <pthread.h>
#include "smiol_types.h"

int SMIOL_async_init(struct SMIOL_context *context);

int SMIOL_async_finalize(struct SMIOL_context *context);

void SMIOL_async_queue_add(struct SMIOL_file *file,
                           struct SMIOL_async_buffer *b);

int SMIOL_async_queue_empty(struct SMIOL_file *file);

struct SMIOL_async_buffer *SMIOL_async_queue_remove(struct SMIOL_file *file);

void SMIOL_async_ticket_lock(struct SMIOL_file *file);

void SMIOL_async_ticket_unlock(struct SMIOL_file *file);

void SMIOL_async_launch_thread(pthread_t **thread,
                               void *(*func)(void *), void *arg);

void SMIOL_async_join_thread(pthread_t **thread);

#endif
