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
	MPI_Fint fcomm; /* Fortran handle to MPI communicator */
	int comm_size;  /* Size of MPI communicator */
	int comm_rank;  /* Rank within MPI communicator */

	int lib_ierr;   /* Library-specific error code */
	int lib_type;   /* From which library the error code originated */
};

struct SMIOL_file {
	struct SMIOL_context *context; /* Context for this file */
	SMIOL_Offset frame; /* Current frame of the file */
#ifdef SMIOL_PNETCDF
	int state; /* parallel-netCDF file state (i.e. Define or data mode) */
	int ncidp; /* parallel-netCDF file handle */
#endif
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
};


/*
 * Return error codes
 */
#include "smiol_codes.inc"

#endif
