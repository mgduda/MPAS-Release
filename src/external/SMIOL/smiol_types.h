/*******************************************************************************
 * SMIOL -- The Simple MPAS I/O Library
 *******************************************************************************/
#ifndef SMIOL_TYPES_H
#define SMIOL_TYPES_H

#include <stdint.h>
#include "mpi.h"


/* If SMIOL_Offset is redefined, interoperable Fortran types and interfaces must also be updated */
typedef int64_t SMIOL_Offset;


#define TRIPLET_SIZE ((size_t)3)


/*
 * Types
 */
struct SMIOL_context {
	MPI_Fint fcomm;   /* Fortran handle to MPI communicator */
	int comm_size;    /* Size of MPI communicator */
	int comm_rank;    /* Rank within MPI communicator */

	int num_io_tasks; /* The number of I/O tasks */
	int io_stride;    /* The stride between I/O tasks in the communicator */

	int lib_ierr;     /* Library-specific error code */
	int lib_type;     /* From which library the error code originated */

	/*
	 * Asynchronous output
	 */
	int async_io_comm;    /* Comm for ranks performing I/O */
	int async_group_comm; /* Comm for group of ranks associated with an I/O rank */
	pthread_mutex_t *mutex;
	pthread_cond_t *cond;

	/*
	 * Checksum for verifying validity of contents of a SMIOL_context struct
	 */
	int checksum;
};

struct SMIOL_file {
	struct SMIOL_context *context; /* Context for this file */
	SMIOL_Offset frame; /* Current frame of the file */
#ifdef SMIOL_PNETCDF
	int state; /* parallel-netCDF file state (i.e. Define or data mode) */
	int ncidp; /* parallel-netCDF file handle */
	int io_task; /* 1 = this task performs I/O calls; 0 = no I/O calls on this task */
	int io_file_comm;
	int io_group_comm;
	int n_reqs;
	int *reqs;
#endif

	/*
	 * Asynchronous output
	 */
	int mode;
	int active;
	pthread_t *writer;
	pthread_mutex_t *mutex;
	pthread_cond_t *cond;
	unsigned long queue_head;
	unsigned long queue_tail;
	struct SMIOL_async_queue *queue;

	/*
	 * Checksum for verifying validity of contents of a SMIOL_file struct
	 */
	int checksum;
};

struct SMIOL_decomp {
	/*
	 * The lists below are structured as follows:
	 *   list[0] - the number of neighbors for which a task sends/recvs
	 *                                                                             |
	 *   list[n] - neighbor task ID                                                | repeated for
	 *   list[n+1] - number of elements, m, to send/recv to/from the neighbor      | each neighbor
	 *   list[n+2 .. n+2+m] - local element IDs to send/recv to/from the neighbor  |
	 *                                                                             |
	 */
	SMIOL_Offset *comp_list;   /* Elements to be sent/received from/on a compute task */
	SMIOL_Offset *io_list;     /* Elements to be sent/received from/on an I/O task */

	struct SMIOL_context *context; /* Context for this decomp */

	size_t io_start;  /* The starting offset on disk for I/O by a task */
	size_t io_count;  /* The number of elements for I/O by a task */

#ifdef SMIOL_AGGREGATION
	MPI_Fint agg_comm;
	size_t n_compute;
	size_t n_compute_agg;
	int *counts;
	int *displs;
#endif
};

struct SMIOL_async_buffer {
	int ierr;
	void *buf;
	size_t bufsize;
#ifdef SMIOL_PNETCDF
	int ncidp;
	int varidp;
	MPI_Offset *mpi_start;
	MPI_Offset *mpi_count;
#endif
	struct SMIOL_async_buffer *next;
};

struct SMIOL_async_queue {
	struct SMIOL_async_buffer *head;
	struct SMIOL_async_buffer *tail;
};

#define SMIOL_ASYNC_QUEUE_INITIALIZER (struct SMIOL_async_queue){NULL, NULL}

/*
 * Return error codes
 */
#include "smiol_codes.inc"

#endif
