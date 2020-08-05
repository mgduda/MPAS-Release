#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "smiol.h"
#include "smiol_utils.h"

#ifdef SMIOL_PNETCDF
#include "pnetcdf.h"
#define PNETCDF_DEFINE_MODE 0
#define PNETCDF_DATA_MODE 1
#endif

#define START_COUNT_READ 0
#define START_COUNT_WRITE 1

/*
 * Local functions
 */
int build_start_count(struct SMIOL_file *file, const char *varname,
                      const struct SMIOL_decomp *decomp,
                      int write_or_read, size_t *element_size, int *ndims,
                      size_t **start, size_t **count);


/********************************************************************************
 *
 * SMIOL_fortran_init
 *
 * Initialize a SMIOL context from Fortran.
 *
 * This function is a simply a wrapper for the SMOIL_init routine that is intended
 * to be called from Fortran. Accordingly, the first argument is of type MPI_Fint
 * (a Fortran integer) rather than MPI_Comm.
 *
 ********************************************************************************/
int SMIOL_fortran_init(MPI_Fint comm, struct SMIOL_context **context)
{
	return SMIOL_init(MPI_Comm_f2c(comm), context);
}


/********************************************************************************
 *
 * SMIOL_init
 *
 * Initialize a SMIOL context.
 *
 * Initializes a SMIOL context, within which decompositions may be defined and
 * files may be read and written. At present, the only input argument is an MPI
 * communicator.
 *
 * Upon successful return the context argument points to a valid SMIOL context;
 * otherwise, it is NULL and an error code other than MPI_SUCCESS is returned.
 *
 * Note: It is assumed that MPI_Init has been called prior to this routine, so
 *       that any use of the provided MPI communicator will be valid.
 *
 ********************************************************************************/
int SMIOL_init(MPI_Comm comm, struct SMIOL_context **context)
{
	MPI_Comm smiol_comm;

	/*
	 * Before dereferencing context below, ensure that the pointer
	 * the context pointer is not NULL
	 */
	if (context == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * We cannot check for every possible invalid argument for comm, but
	 * at least we can verify that the communicator is not MPI_COMM_NULL
	 */
	if (comm == MPI_COMM_NULL) {
		/* Nullifying (*context) here may result in a memory leak, but this
		 * seems better than disobeying the stated behavior of returning
		 * a NULL context upon failure
		 */
		(*context) = NULL;

		return SMIOL_INVALID_ARGUMENT;
	}

	*context = (struct SMIOL_context *)malloc(sizeof(struct SMIOL_context));
	if ((*context) == NULL) {
		return SMIOL_MALLOC_FAILURE;
	}

	/*
	 * Initialize context
	 */
	(*context)->lib_ierr = 0;
	(*context)->lib_type = SMIOL_LIBRARY_UNKNOWN;

	/*
	 * Make a duplicate of the MPI communicator for use by SMIOL
	 */
	if (MPI_Comm_dup(comm, &smiol_comm) != MPI_SUCCESS) {
		free((*context));
		(*context) = NULL;
		return SMIOL_MPI_ERROR;
	}
	(*context)->fcomm = MPI_Comm_c2f(smiol_comm);

	if (MPI_Comm_size(smiol_comm, &((*context)->comm_size)) != MPI_SUCCESS) {
		free((*context));
		(*context) = NULL;
		return SMIOL_MPI_ERROR;
	}

	if (MPI_Comm_rank(smiol_comm, &((*context)->comm_rank)) != MPI_SUCCESS) {
		free((*context));
		(*context) = NULL;
		return SMIOL_MPI_ERROR;
	}

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_finalize
 *
 * Finalize a SMIOL context.
 *
 * Finalizes a SMIOL context and frees all memory in the SMIOL_context instance.
 * After this routine is called, no other SMIOL routines that make reference to
 * the finalized context should be called.
 *
 ********************************************************************************/
int SMIOL_finalize(struct SMIOL_context **context)
{
	MPI_Comm smiol_comm;

	/*
	 * If the pointer to the context pointer is NULL, assume we have nothing
	 * to do and declare success
	 */
	if (context == NULL) {
		return SMIOL_SUCCESS;
	}

	if ((*context) == NULL) {
		return SMIOL_SUCCESS;
	}

	smiol_comm = MPI_Comm_f2c((*context)->fcomm);
	if (MPI_Comm_free(&smiol_comm) != MPI_SUCCESS) {
		free((*context));
		(*context) = NULL;
		return SMIOL_MPI_ERROR;
	}

	free((*context));
	(*context) = NULL;

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_inquire
 *
 * Inquire about a SMIOL context.
 *
 * Detailed description.
 *
 ********************************************************************************/
int SMIOL_inquire(void)
{
	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_open_file
 *
 * Opens a file within a SMIOL context.
 *
 * Depending on the specified file mode, creates or opens the file specified
 * by filename within the provided SMIOL context.
 *
 * Upon successful completion, SMIOL_SUCCESS is returned, and the file handle
 * argument will point to a valid file handle and the current frame for the
 * file will be set to zero. Otherwise, the file handle is NULL and an error
 * code other than SMIOL_SUCCESS is returned.
 *
 ********************************************************************************/
int SMIOL_open_file(struct SMIOL_context *context, const char *filename, int mode, struct SMIOL_file **file)
{
#ifdef SMIOL_PNETCDF
	int ierr;
#endif

	/*
	 * Before dereferencing file below, ensure that the pointer
	 * the file pointer is not NULL
	 */
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Check that context is valid
	 */
	if (context == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	*file = (struct SMIOL_file *)malloc(sizeof(struct SMIOL_file));
	if ((*file) == NULL) {
		return SMIOL_MALLOC_FAILURE;
	}

	/*
	 * Save pointer to context for this file
	 */
	(*file)->context = context;
	(*file)->frame = (SMIOL_Offset) 0;

	if (mode & SMIOL_FILE_CREATE) {
#ifdef SMIOL_PNETCDF
		if ((ierr = ncmpi_create(MPI_Comm_f2c(context->fcomm), filename,
					(NC_64BIT_DATA | NC_CLOBBER), MPI_INFO_NULL,
					&((*file)->ncidp))) != NC_NOERR) {
			free((*file));
			(*file) = NULL;
			context->lib_type = SMIOL_LIBRARY_PNETCDF;
			context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		} else {
			(*file)->state = PNETCDF_DEFINE_MODE;
		}
#endif
	}
	else if (mode & SMIOL_FILE_WRITE) {
#ifdef SMIOL_PNETCDF
		if ((ierr = ncmpi_open(MPI_Comm_f2c(context->fcomm), filename,
				NC_WRITE, MPI_INFO_NULL, &((*file)->ncidp))) != NC_NOERR) {
			free((*file));
			(*file) = NULL;
			context->lib_type = SMIOL_LIBRARY_PNETCDF;
			context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		} else {
			(*file)->state = PNETCDF_DATA_MODE;
		}
#endif
	}
	else if (mode & SMIOL_FILE_READ) {
#ifdef SMIOL_PNETCDF
		if ((ierr = ncmpi_open(MPI_Comm_f2c(context->fcomm), filename,
				NC_NOWRITE, MPI_INFO_NULL, &((*file)->ncidp))) != NC_NOERR) {
			free((*file));
			(*file) = NULL;
			context->lib_type = SMIOL_LIBRARY_PNETCDF;
			context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		} else {
			(*file)->state = PNETCDF_DATA_MODE;
		}
#endif
	}
	else {
		free((*file));
		(*file) = NULL;
		return SMIOL_INVALID_ARGUMENT;
	}

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_close_file
 *
 * Closes a file within a SMIOL context.
 *
 * Closes the file associated with the provided file handle. Upon successful
 * completion, SMIOL_SUCCESS is returned, the file will be closed, and all memory
 * that is uniquely associated with the file handle will be deallocated.
 * Otherwise, an error code other than SMIOL_SUCCESS will be returned.
 *
 ********************************************************************************/
int SMIOL_close_file(struct SMIOL_file **file)
{
#ifdef SMIOL_PNETCDF
	int ierr;
#endif

	/*
	 * If the pointer to the file pointer is NULL, assume we have nothing
	 * to do and declare success
	 */
	if (file == NULL) {
		return SMIOL_SUCCESS;
	}

#ifdef SMIOL_PNETCDF
	if ((ierr = ncmpi_close((*file)->ncidp)) != NC_NOERR) {
		((*file)->context)->lib_type = SMIOL_LIBRARY_PNETCDF;
		((*file)->context)->lib_ierr = ierr;
		free((*file));
		(*file) = NULL;
		return SMIOL_LIBRARY_ERROR;
	}
#endif

	free((*file));
	(*file) = NULL;

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_define_dim
 *
 * Defines a new dimension in a file.
 *
 * Defines a dimension with the specified name and size in the file associated
 * with the file handle. If a negative value is provided for the size argument,
 * the dimension will be defined as an unlimited or record dimension.
 *
 * Upon successful completion, SMIOL_SUCCESS is returned; otherwise, an error
 * code is returned.
 *
 ********************************************************************************/
int SMIOL_define_dim(struct SMIOL_file *file, const char *dimname, SMIOL_Offset dimsize)
{
#ifdef SMIOL_PNETCDF
	int dimidp;
	int ierr;
	MPI_Offset len;
#endif

	/*
	 * Check that file handle is valid
	 */
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Check that dimension name is valid
	 */
	if (dimname == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

#ifdef SMIOL_PNETCDF
	/*
	 * The parallel-netCDF library does not permit zero-length dimensions
	 */
	if (dimsize == (SMIOL_Offset)0) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Handle unlimited / record dimension specifications
	 */
	if (dimsize < (SMIOL_Offset)0) {
		len = NC_UNLIMITED;
	}
	else {
		len = (MPI_Offset)dimsize;
	}

	/*
	 * If the file is in data mode, then switch it to define mode
	 */
	if (file->state == PNETCDF_DATA_MODE) {
		if ((ierr = ncmpi_redef(file->ncidp)) != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
		file->state = PNETCDF_DEFINE_MODE;
	}

	if ((ierr = ncmpi_def_dim(file->ncidp, dimname, len, &dimidp)) != NC_NOERR) {
		file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
		file->context->lib_ierr = ierr;
		return SMIOL_LIBRARY_ERROR;
	}
#endif

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_inquire_dim
 *
 * Inquires about an existing dimension in a file.
 *
 * Inquire about the size of an existing dimension and optionally inquire if the
 * given dimension is the unlimited dimension or not. If dimsize is a non-NULL
 * pointer then the dimension size will be returned in dimsize. For unlimited
 * dimensions, the current size of the dimension is returned; future writes of
 * additional records to a file can lead to different return sizes for
 * unlimited dimensions.
 *
 * If is_unlimited is a non-NULL pointer and if the inquired dimension is the
 * unlimited dimension, is_unlimited will be set to 1; if the inquired
 * dimension is not the unlimited dimension then is_unlimited will be set to 0.
 *
 * Upon successful completion, SMIOL_SUCCESS is returned; otherwise, an error
 * code is returned.
 *
 ********************************************************************************/
int SMIOL_inquire_dim(struct SMIOL_file *file, const char *dimname,
                      SMIOL_Offset *dimsize, int *is_unlimited)
{
#ifdef SMIOL_PNETCDF
	int dimidp;
	int ierr;
	MPI_Offset len;
#endif
	/*
	 * Check that file handle is valid
	 */
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Check that dimension name is valid
	 */
	if (dimname == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Check that dimension size is not NULL
	 */
	if (dimsize == NULL && is_unlimited == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	if (dimsize != NULL) {
		(*dimsize) = (SMIOL_Offset)0;   /* Default dimension size if no library provides a value */
	}

	if (is_unlimited != NULL) {
		(*is_unlimited) = 0; /* Return 0 if no library provides a value */
	}

#ifdef SMIOL_PNETCDF
	if ((ierr = ncmpi_inq_dimid(file->ncidp, dimname, &dimidp)) != NC_NOERR) {
		(*dimsize) = (SMIOL_Offset)(-1);  /* TODO: should there be a well-defined invalid size? */
		file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
		file->context->lib_ierr = ierr;
		return SMIOL_LIBRARY_ERROR;
	}

	/*
	 * Inquire about dimsize
	 */
	if (dimsize != NULL) {
		ierr = ncmpi_inq_dimlen(file->ncidp, dimidp, &len);
		if (ierr != NC_NOERR) {
			(*dimsize) = (SMIOL_Offset)(-1);  /* TODO: should there be a well-defined invalid size? */
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}

		(*dimsize) = (SMIOL_Offset)len;
	}


	/*
	 * Inquire if this dimension is the unlimited dimension
	 */
	if (is_unlimited != NULL) {
		int unlimdimidp;
		ierr = ncmpi_inq_unlimdim(file->ncidp, &unlimdimidp);
		if (ierr != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
		if (unlimdimidp == dimidp) {
			(*is_unlimited) = 1;
		} else {
			(*is_unlimited) = 0; // Not the unlimited dim
		}
	}
#endif

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_define_var
 *
 * Defines a new variable in a file.
 *
 * Defines a variable with the specified name, type, and dimensions in an open
 * file pointed to by the file argument. The varname and dimnames arguments
 * are expected to be null-terminated strings, except if the variable has zero
 * dimensions, in which case the dimnames argument may be a NULL pointer.
 *
 * Upon successful completion, SMIOL_SUCCESS is returned; otherwise, an error
 * code is returned.
 *
 ********************************************************************************/
int SMIOL_define_var(struct SMIOL_file *file, const char *varname, int vartype, int ndims, const char **dimnames)
{
#ifdef SMIOL_PNETCDF
	int *dimids;
	int ierr;
	int i;
	nc_type xtype;
	int varidp;
#endif

	/*
	 * Check that file handle is valid
	 */
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Check that variable name is valid
	 */
	if (varname == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Check that the variable type is valid - handled below in a library-specific way...
	 */

	/*
	 * Check that variable dimension names are valid
	 */
	if (dimnames == NULL && ndims > 0) {
		return SMIOL_INVALID_ARGUMENT;
	}

#ifdef SMIOL_PNETCDF
	dimids = (int *)malloc(sizeof(int) * (size_t)ndims);
	if (dimids == NULL) {
		return SMIOL_MALLOC_FAILURE;
	}

	/*
	 * Build a list of dimension IDs
	 */
	for (i=0; i<ndims; i++) {
		if ((ierr = ncmpi_inq_dimid(file->ncidp, dimnames[i], &dimids[i])) != NC_NOERR) {
			free(dimids);
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
	}

	/*
	 * Translate SMIOL variable type to parallel-netcdf type
	 */
	switch (vartype) {
		case SMIOL_REAL32:
			xtype = NC_FLOAT;
			break;
		case SMIOL_REAL64:
			xtype = NC_DOUBLE;
			break;
		case SMIOL_INT32:
			xtype = NC_INT;
			break;
		case SMIOL_CHAR:
			xtype = NC_CHAR;
			break;
		default:
			free(dimids);
			return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * If the file is in data mode, then switch it to define mode
	 */
	if (file->state == PNETCDF_DATA_MODE) {
		if ((ierr = ncmpi_redef(file->ncidp)) != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
		file->state = PNETCDF_DEFINE_MODE;
	}

	/*
	 * Define the variable
	 */
	if ((ierr = ncmpi_def_var(file->ncidp, varname, xtype, ndims, dimids, &varidp)) != NC_NOERR) {
		free(dimids);
		file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
		file->context->lib_ierr = ierr;
		return SMIOL_LIBRARY_ERROR;
	}

	free(dimids);
#endif

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_inquire_var
 *
 * Inquires about an existing variable in a file.
 *
 * Inquires about a variable in a file, and optionally returns the type
 * of the variable, the dimensionality of the variable, and the names of
 * the dimensions of the variable. Which properties of the variable to return
 * (type, dimensionality, or dimension names) is indicated by the status of
 * the pointers for the corresponding properties: if the pointer is a non-NULL
 * pointer, the property will be set upon successful completion of this routine.
 *
 * If the names of a variable's dimensions are requested (by providing a non-NULL
 * actual argument for dimnames), the size of the dimnames array must be at least
 * the number of dimensions in the variable, and each character string pointed
 * to by an element of dimnames must be large enough to accommodate the corresponding
 * dimension name.
 *
 ********************************************************************************/
int SMIOL_inquire_var(struct SMIOL_file *file, const char *varname, int *vartype, int *ndims, char **dimnames)
{
#ifdef SMIOL_PNETCDF
	int *dimids;
	int varidp;
	int ierr;
	int i;
	int xtypep;
	int ndimsp;
#endif

	/*
	 * Check that file handle is valid
	 */
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Check that variable name is valid
	 */
	if (varname == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * If all output arguments are NULL, we can return early
	 */
	if (vartype == NULL && ndims == NULL && dimnames == NULL) {
		return SMIOL_SUCCESS;
	}

	/*
	 * Provide default values for output arguments in case
	 * no library-specific below is active
	 */
	if (vartype != NULL) {
		*vartype = SMIOL_UNKNOWN_VAR_TYPE;
	}
	if (ndims != NULL) {
		*ndims = 0;
	}

#ifdef SMIOL_PNETCDF
	/*
	 * Get variable ID
	 */
	if ((ierr = ncmpi_inq_varid(file->ncidp, varname, &varidp)) != NC_NOERR) {
		file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
		file->context->lib_ierr = ierr;
		return SMIOL_LIBRARY_ERROR;
	}

	/*
	 * If requested, inquire about variable type
	 */
	if (vartype != NULL) {
		if ((ierr = ncmpi_inq_vartype(file->ncidp, varidp, &xtypep)) != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}

		/* Convert parallel-netCDF variable type to SMIOL variable type */
		switch (xtypep) {
			case NC_FLOAT:
				*vartype = SMIOL_REAL32;
				break;
			case NC_DOUBLE:
				*vartype = SMIOL_REAL64;
				break;
			case NC_INT:
				*vartype = SMIOL_INT32;
				break;
			case NC_CHAR:
				*vartype = SMIOL_CHAR;
				break;
			default:
				*vartype = SMIOL_UNKNOWN_VAR_TYPE;
		}
	}

	/*
	 * All remaining properties will require the number of dimensions
	 */
	if (ndims != NULL || dimnames != NULL) {
		if ((ierr = ncmpi_inq_varndims(file->ncidp, varidp, &ndimsp)) != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
	}

	/*
	 * If requested, inquire about dimensionality
	 */
	if (ndims != NULL) {
		*ndims = ndimsp;
	}

	/*
	 * If requested, inquire about dimension names
	 */
	if (dimnames != NULL) {
		dimids = (int *)malloc(sizeof(int) * (size_t)ndimsp);
		if (dimids == NULL) {
			return SMIOL_MALLOC_FAILURE;
		}

		if ((ierr = ncmpi_inq_vardimid(file->ncidp, varidp, dimids)) != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			free(dimids);
			return SMIOL_LIBRARY_ERROR;
		}

		for (i = 0; i < ndimsp; i++) {
			if (dimnames[i] == NULL) {
				return SMIOL_INVALID_ARGUMENT;
			}
			if ((ierr = ncmpi_inq_dimname(file->ncidp, dimids[i], dimnames[i])) != NC_NOERR) {
				file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
				file->context->lib_ierr = ierr;
				free(dimids);
				return SMIOL_LIBRARY_ERROR;
			}
		}

		free(dimids);
	}
#endif

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_put_var
 *
 * Writes a variable to a file.
 *
 * Given a pointer to a SMIOL file that was previously opened with write access
 * and the name of a variable previously defined in the file with a call to
 * SMIOL_define_var, this routine will write the contents of buf to the variable
 * according to the decomposition described by decomp.
 *
 * If decomp is not NULL, the variable is assumed to be decomposed across MPI
 * ranks, and all ranks with non-zero-sized partitions of the variable must
 * provide a valid buffer. For decomposed variables, all MPI ranks must provide
 * a non-NULL decomp, regardless of whether a rank has a non-zero-sized
 * partition of the variable.
 *
 * If the variable is not decomposed -- that is, all ranks store identical
 * values for the entire variable -- all MPI ranks must provide a NULL pointer
 * for the decomp argument. As currently implemented, this routine will write
 * the buffer for MPI rank 0 to the variable; however, this behavior should not
 * be relied on.
 *
 * If the variable has been successfully written to the file, SMIOL_SUCCESS will
 * be returned. Otherwise, an error code indicating the nature of the failure
 * will be returned.
 *
 ********************************************************************************/
int SMIOL_put_var(struct SMIOL_file *file, const char *varname,
                  const struct SMIOL_decomp *decomp, const void *buf)
{
	int ierr;
	int ndims;
	size_t element_size;
	void *out_buf = NULL;
	size_t *start;
	size_t *count;

	/*
	 * Basic checks on arguments
	 */
	if (file == NULL || varname == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Work out the start[] and count[] arrays for writing this variable
	 * in parallel
	 */
	ierr = build_start_count(file, varname, decomp,
	                         START_COUNT_WRITE, &element_size, &ndims,
	                         &start, &count);
	if (ierr != SMIOL_SUCCESS) {
		return ierr;
	}

	/*
	 * Communicate elements of this field from MPI ranks that compute those
	 * elements to MPI ranks that write those elements. This only needs to
	 * be done for decomposed variables.
	 */
	if (decomp) {
		out_buf = malloc(element_size * decomp->io_count);
		if (out_buf == NULL) {
			free(start);
			free(count);

			return SMIOL_MALLOC_FAILURE;
		}

		ierr = transfer_field(decomp, SMIOL_COMP_TO_IO,
		                      element_size, buf, out_buf);
		if (ierr != SMIOL_SUCCESS) {
			free(start);
			free(count);
			free(out_buf);
			return ierr;
		}
	}

	/*
	 * Write out_buf
	 */
#ifdef SMIOL_PNETCDF
	{
		int j;
		int varidp;
		const void *buf_p;
		MPI_Offset *mpi_start;
		MPI_Offset *mpi_count;

		if (file->state == PNETCDF_DEFINE_MODE) {
			if ((ierr = ncmpi_enddef(file->ncidp)) != NC_NOERR) {
				file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
				file->context->lib_ierr = ierr;

				if (decomp) {
					free(out_buf);
				}
				free(start);
				free(count);

				return SMIOL_LIBRARY_ERROR;
			}
			file->state = PNETCDF_DATA_MODE;
		}

		ierr = ncmpi_inq_varid(file->ncidp, varname, &varidp);
		if (ierr != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;

			if (decomp) {
				free(out_buf);
			}
			free(start);
			free(count);

			return SMIOL_LIBRARY_ERROR;
		}

		if (decomp) {
			buf_p = out_buf;
		} else {
			buf_p = buf;
		}

		mpi_start = malloc(sizeof(MPI_Offset) * (size_t)ndims);
		if (mpi_start == NULL) {
			free(start);
			free(count);

			return SMIOL_MALLOC_FAILURE;
		}

		mpi_count = malloc(sizeof(MPI_Offset) * (size_t)ndims);
		if (mpi_count == NULL) {
			free(start);
			free(count);
			free(mpi_start);

			return SMIOL_MALLOC_FAILURE;
		}

		for (j = 0; j < ndims; j++) {
			mpi_start[j] = (MPI_Offset)start[j];
			mpi_count[j] = (MPI_Offset)count[j];
		}

		ierr = ncmpi_put_vara_all(file->ncidp,
		                          varidp,
		                          mpi_start, mpi_count,
		                          buf_p,
		                          0, MPI_DATATYPE_NULL);

		free(mpi_start);
		free(mpi_count);

		if (ierr != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;

			if (decomp) {
				free(out_buf);
			}
			free(start);
			free(count);

			return SMIOL_LIBRARY_ERROR;
		}
	}
#endif

	/*
	 * Free up memory before returning
	 */
	if (decomp) {
		free(out_buf);
	}

	free(start);
	free(count);

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_get_var
 *
 * Reads a variable from a file.
 *
 * Given a pointer to a SMIOL file and the name of a variable previously defined
 * in the file, this routine will read the contents of the variable into buf
 * according to the decomposition described by decomp.
 *
 * If decomp is not NULL, the variable is assumed to be decomposed across MPI
 * ranks, and all ranks with non-zero-sized partitions of the variable must
 * provide a valid buffer. For decomposed variables, all MPI ranks must provide
 * a non-NULL decomp, regardless of whether a rank has a non-zero-sized
 * partition of the variable.
 *
 * If the variable is not decomposed -- that is, all ranks load identical
 * values for the entire variable -- all MPI ranks must provide a NULL pointer
 * for the decomp argument.
 *
 * If the variable has been successfully read from the file, SMIOL_SUCCESS will
 * be returned. Otherwise, an error code indicating the nature of the failure
 * will be returned.
 *
 ********************************************************************************/
int SMIOL_get_var(struct SMIOL_file *file, const char *varname,
                  const struct SMIOL_decomp *decomp, void *buf)
{
	int ierr;
	int ndims;
	size_t element_size;
	void *in_buf = NULL;
	size_t *start;
	size_t *count;

	/*
	 * Basic checks on arguments
	 */
	if (file == NULL || varname == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Work out the start[] and count[] arrays for reading this variable
	 * in parallel
	 */
	ierr = build_start_count(file, varname, decomp,
	                         START_COUNT_READ, &element_size, &ndims,
	                         &start, &count);
	if (ierr != SMIOL_SUCCESS) {
		return ierr;
	}

	/*
	 * If this variable is decomposed, allocate a buffer into which
	 * the variable will be read using the I/O decomposition; later,
	 * elements this buffer will be transferred to MPI ranks that compute
	 * on those elements
	 */
	if (decomp) {
		in_buf = malloc(element_size * decomp->io_count);
		if (in_buf == NULL) {
			free(start);
			free(count);

			return SMIOL_MALLOC_FAILURE;
		}

#ifndef SMIOL_PNETCDF
		/*
		 * If no file library provides values for the memory pointed to
		 * by in_buf, the transfer_field call later will transfer
		 * garbage to the output buffer; to avoid returning
		 * non-deterministic values to the caller in this case,
		 * initialize in_buf.
		 */
		memset(in_buf, 0, element_size * decomp->io_count);
		
#endif
	}

	/*
	 * Read in_buf
	 */
#ifdef SMIOL_PNETCDF
	{
		int j;
		int varidp;
		void *buf_p;
		MPI_Offset *mpi_start;
		MPI_Offset *mpi_count;

		if (file->state == PNETCDF_DEFINE_MODE) {
			if ((ierr = ncmpi_enddef(file->ncidp)) != NC_NOERR) {
				file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
				file->context->lib_ierr = ierr;

				if (decomp) {
					free(in_buf);
				}
				free(start);
				free(count);

				return SMIOL_LIBRARY_ERROR;
			}
			file->state = PNETCDF_DATA_MODE;
		}

		ierr = ncmpi_inq_varid(file->ncidp, varname, &varidp);
		if (ierr != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;

			if (decomp) {
				free(in_buf);
			}
			free(start);
			free(count);

			return SMIOL_LIBRARY_ERROR;
		}

		if (decomp) {
			buf_p = in_buf;
		} else {
			buf_p = buf;
		}

		mpi_start = malloc(sizeof(MPI_Offset) * (size_t)ndims);
		if (mpi_start == NULL) {
			free(start);
			free(count);

			return SMIOL_MALLOC_FAILURE;
		}

		mpi_count = malloc(sizeof(MPI_Offset) * (size_t)ndims);
		if (mpi_count == NULL) {
			free(start);
			free(count);
			free(mpi_start);

			return SMIOL_MALLOC_FAILURE;
		}

		for (j = 0; j < ndims; j++) {
			mpi_start[j] = (MPI_Offset)start[j];
			mpi_count[j] = (MPI_Offset)count[j];
		}

		ierr = ncmpi_get_vara_all(file->ncidp,
		                          varidp,
		                          mpi_start, mpi_count,
		                          buf_p,
		                          0, MPI_DATATYPE_NULL);

		free(mpi_start);
		free(mpi_count);

		if (ierr != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;

			if (decomp) {
				free(in_buf);
			}
			free(start);
			free(count);

			return SMIOL_LIBRARY_ERROR;
		}
	}
#endif

	/*
	 * Free start/count arrays
	 */
	free(start);
	free(count);

	/*
	 * Communicate elements of this field from MPI ranks that read those
	 * elements to MPI ranks that compute those elements. This only needs to
	 * be done for decomposed variables.
	 */
	if (decomp) {
		ierr = transfer_field(decomp, SMIOL_IO_TO_COMP,
		                      element_size, in_buf, buf);
		free(in_buf);

		if (ierr != SMIOL_SUCCESS) {
			return ierr;
		}
	}

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_define_att
 *
 * Defines a new attribute in a file.
 *
 * Defines a new attribute for a variable if varname is not NULL,
 * or a global attribute otherwise. The type of the attribute must be one
 * of SMIOL_REAL32, SMIOL_REAL64, SMIOL_INT32, or SMIOL_CHAR.
 *
 * If the attribute has been successfully defined for the variable or file,
 * SMIOL_SUCCESS is returned.
 *
 ********************************************************************************/
int SMIOL_define_att(struct SMIOL_file *file, const char *varname,
                     const char *att_name, int att_type, const void *att)
{
#ifdef SMIOL_PNETCDF
	int ierr;
	int varidp;
	nc_type xtype;
#endif

	/*
	 * Check validity of arguments
	 */
	if (file == NULL || att_name == NULL || att == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Checks for valid attribute type are handled in library-specific
	 * code, below
	 */

#ifdef SMIOL_PNETCDF
	/*
	 * If varname was provided, get the variable ID; else, the attribute
	 * is a global attribute not associated with a specific variable
	 */
	if (varname != NULL) {
		ierr = ncmpi_inq_varid(file->ncidp, varname, &varidp);
		if (ierr != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
	} else {
		varidp = NC_GLOBAL;
	}

	/*
	 * Translate SMIOL variable type to parallel-netcdf type
	 */
	switch (att_type) {
		case SMIOL_REAL32:
			xtype = NC_FLOAT;
			break;
		case SMIOL_REAL64:
			xtype = NC_DOUBLE;
			break;
		case SMIOL_INT32:
			xtype = NC_INT;
			break;
		case SMIOL_CHAR:
			xtype = NC_CHAR;
			break;
		default:
			return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * If the file is in data mode, then switch it to define mode
	 */
	if (file->state == PNETCDF_DATA_MODE) {
		if ((ierr = ncmpi_redef(file->ncidp)) != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
		file->state = PNETCDF_DEFINE_MODE;
	}

	/*
	 * Add the attribute to the file
	 */
	if (att_type == SMIOL_CHAR) {
		ierr = ncmpi_put_att(file->ncidp, varidp, att_name, xtype,
		                     (MPI_Offset)strlen(att), (const char *)att);
	} else {
		ierr = ncmpi_put_att(file->ncidp, varidp, att_name, xtype,
		                     (MPI_Offset)1, (const char *)att);
	}
	if (ierr != NC_NOERR) {
		file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
		file->context->lib_ierr = ierr;
		return SMIOL_LIBRARY_ERROR;
	}
#endif

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_inquire_att
 *
 * Inquires about an attribute in a file.
 *
 * Inquires about a variable attribute if varname is not NULL, or a global
 * attribute otherwise.
 *
 * If the requested attribute is found, SMIOL_SUCCESS is returned and the memory
 * pointed to by the att argument will contain the attribute value.
 *
 * For character string attributes, no bytes beyond the length of the attribute
 * in the file will be modified in the att argument, and no '\0' character will
 * be added. Therefore, calling code may benefit from initializing character
 * strings before calling this routine.
 *
 * If SMIOL was not compiled with support for any file library, the att_type
 * output argument will always be set to SMIOL_UNKNOWN_VAR_TYPE, and the att_len
 * output argument will always be set to -1; the value of the att output
 * argument will be unchanged.
 *
 ********************************************************************************/
int SMIOL_inquire_att(struct SMIOL_file *file, const char *varname,                                                                  
                      const char *att_name, int *att_type,
                      SMIOL_Offset *att_len, void *att)
{
#ifdef SMIOL_PNETCDF
	int ierr;
	int varidp;
	nc_type xtypep;
	MPI_Offset lenp;
#endif

	/*
	 * Check validity of arguments
	 */
	if (file == NULL || att_name == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	/*
	 * Set output arguments in case no library sets them later
	 */
	if (att_len != NULL) {
		*att_len = (SMIOL_Offset)-1;
	}

	if (att_type != NULL) {
		*att_type = SMIOL_UNKNOWN_VAR_TYPE;
	}

#ifdef SMIOL_PNETCDF
	/*
	 * If varname was provided, get the variable ID; else, the inquiry is
	 * is for a global attribute not associated with a specific variable
	 */
	if (varname != NULL) {
		ierr = ncmpi_inq_varid(file->ncidp, varname, &varidp);
		if (ierr != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
	} else {
		varidp = NC_GLOBAL;
	}

	/*
	 * Inquire about attribute type and length
	 */
	if (att_type != NULL || att_len != NULL) {
		ierr = ncmpi_inq_att(file->ncidp, varidp, att_name,
		                     &xtypep, &lenp);
		if (ierr != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}

		if (att_type != NULL) {
			/* Convert parallel-netCDF type to SMIOL type */
			switch (xtypep) {
				case NC_FLOAT:
					*att_type = SMIOL_REAL32;
					break;
				case NC_DOUBLE:
					*att_type = SMIOL_REAL64;
					break;
				case NC_INT:
					*att_type = SMIOL_INT32;
					break;
				case NC_CHAR:
					*att_type = SMIOL_CHAR;
					break;
				default:
					*att_type = SMIOL_UNKNOWN_VAR_TYPE;
			}
		}

		if (att_len != NULL) {
			*att_len = lenp;
		}
	}


	/*
	 * Inquire about attribute value if requested
	 */
	if (att != NULL) {
		ierr = ncmpi_get_att(file->ncidp, varidp, att_name, att);
		if (ierr != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
	}
#endif

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_sync_file
 *
 * Forces all in-memory data to be flushed to disk.
 *
 * Upon success, all in-memory data for the file associatd with the file
 * handle will be flushed to the file system and SMIOL_SUCCESS will be
 * returned; otherwise, an error code is returned.
 *
 ********************************************************************************/
int SMIOL_sync_file(struct SMIOL_file *file)
{
#ifdef SMIOL_PNETCDF
	int ierr;
#endif

	/*
	 * Check that file is valid
	 */
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

#ifdef SMIOL_PNETCDF
	/*
	 * If the file is in define mode then switch it into data mode
	 */
	if (file->state == PNETCDF_DEFINE_MODE) {
		if ((ierr = ncmpi_enddef(file->ncidp)) != NC_NOERR) {
			file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
			file->context->lib_ierr = ierr;
			return SMIOL_LIBRARY_ERROR;
		}
		file->state = PNETCDF_DATA_MODE;
	}

	if ((ierr = ncmpi_sync(file->ncidp)) != NC_NOERR) {
		file->context->lib_type = SMIOL_LIBRARY_PNETCDF;
		file->context->lib_ierr = ierr;
		return SMIOL_LIBRARY_ERROR;
	}
#endif

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * SMIOL_error_string
 *
 * Returns an error string for a specified error code.
 *
 * Returns an error string corresponding to a SMIOL error code. If the error code is
 * SMIOL_LIBRARY_ERROR and a valid SMIOL context is available, the SMIOL_lib_error_string
 * function should be called instead. The error string is null-terminated, but it
 * does not contain a newline character.
 *
 ********************************************************************************/
const char *SMIOL_error_string(int errno)
{
	switch (errno) {
	case SMIOL_SUCCESS:
		return "Success!";
	case SMIOL_MALLOC_FAILURE:
		return "malloc returned a null pointer";
	case SMIOL_INVALID_ARGUMENT:
		return "invalid subroutine argument";
	case SMIOL_MPI_ERROR:
		return "internal MPI call failed";
	case SMIOL_FORTRAN_ERROR:
		return "Fortran wrapper detected an inconsistency in C return values";
	case SMIOL_LIBRARY_ERROR:
		return "bad return code from a library call";
	case SMIOL_WRONG_ARG_TYPE:
		return "argument is of the wrong type";
	case SMIOL_INSUFFICIENT_ARG:
		return "argument is of insufficient size";
	default:
		return "Unknown error";
	}
}


/********************************************************************************
 *
 * SMIOL_lib_error_string
 *
 * Returns an error string for a third-party library called by SMIOL.
 *
 * Returns an error string corresponding to an error that was generated by
 * a third-party library that was called by SMIOL. The library that was the source
 * of the error, as well as the library-specific error code, are retrieved from
 * a SMIOL context. If successive library calls resulted in errors, only the error
 * string for the last of these errors will be returned. The error string is
 * null-terminated, but it does not contain a newline character.
 *
 ********************************************************************************/
const char *SMIOL_lib_error_string(struct SMIOL_context *context)
{
	if (context == NULL) {
		return "SMIOL_context argument is a NULL pointer";
	}

	switch (context->lib_type) {
#ifdef SMIOL_PNETCDF
	case SMIOL_LIBRARY_PNETCDF:
		return ncmpi_strerror(context->lib_ierr);
#endif
	default:
		return "Could not find matching library for the source of the error";
	}
}


/********************************************************************************
 *
 * SMIOL_set_option
 *
 * Sets an option for the SMIOL library.
 *
 * Detailed description.
 *
 ********************************************************************************/
int SMIOL_set_option(void)
{
	return SMIOL_SUCCESS;
}

/********************************************************************************
 *
 * SMIOL_set_frame
 *
 * Set the frame for the unlimited dimension for an open file
 *
 * For an open SMIOL file handle, set the frame for the unlimited dimension.
 * After setting the frame for a file, writing to a variable that is
 * dimensioned by the unlimited dimension will write to the last set frame,
 * overwriting any current data that maybe present in that frame.
 *
 * SMIOL_SUCCESS will be returned if the frame is successfully set otherwise an
 * error will return.
 *
 ********************************************************************************/
int SMIOL_set_frame(struct SMIOL_file *file, SMIOL_Offset frame)
{
	if (file == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}
	file->frame = frame;
	return SMIOL_SUCCESS;
}

/********************************************************************************
 *
 * SMIOL_get_frame
 *
 * Return the current frame of an open file
 *
 * Get the current frame of an open file. Upon success, SMIOL_SUCCESS will be
 * returned, otherwise an error will be returned.
 *
 ********************************************************************************/
int SMIOL_get_frame(struct SMIOL_file *file, SMIOL_Offset *frame)
{
	if (file == NULL || frame == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}
	*frame = file->frame;
	return SMIOL_SUCCESS;
}


/*******************************************************************************
 *
 * SMIOL_create_decomp
 *
 * Creates a mapping between compute elements and I/O elements.
 *
 * Given arrays of global element IDs that each task computes, the number of I/O
 * tasks, and the stride between I/O tasks, this routine works out a mapping of
 * elements between compute and I/O tasks.
 *
 * If all input arguments are determined to be valid and if the routine is
 * successful in working out a mapping, the decomp pointer is allocated and
 * given valid contents, and SMIOL_SUCCESS is returned; otherwise a non-success
 * error code is returned and the decomp pointer is NULL.
 *
 *******************************************************************************/
int SMIOL_create_decomp(struct SMIOL_context *context,
                        size_t n_compute_elements, SMIOL_Offset *compute_elements,
                        int num_io_tasks, int io_stride,
                        struct SMIOL_decomp **decomp)
{
	size_t i;
	size_t n_io_elements, n_io_elements_global;
	size_t io_start, io_count;
	SMIOL_Offset *io_elements;
	MPI_Comm comm;
	MPI_Datatype dtype;
	int ierr;


	/*
	 * Minimal check on the validity of arguments
	 */
	if (context == NULL) {
		return SMIOL_INVALID_ARGUMENT;
	}

	if (compute_elements == NULL && n_compute_elements != 0) {
		return SMIOL_INVALID_ARGUMENT;
	}

	comm = MPI_Comm_f2c(context->fcomm);

	/*
	 * Figure out MPI_Datatype for size_t... there must be a better way...
	 */
	switch (sizeof(size_t)) {
		case sizeof(uint64_t):
			dtype = MPI_UINT64_T;
			break;
		case sizeof(uint32_t):
			dtype = MPI_UINT32_T;
			break;
		case sizeof(uint16_t):
			dtype = MPI_UINT16_T;
			break;
		default:
			return SMIOL_MPI_ERROR;
	}

	/*
	 * Based on the number of compute elements for each task, determine
	 * the total number of elements across all tasks for I/O. The assumption
	 * is that the number of elements to read/write is equal to the size of
	 * the set of compute elements.
	 */
	n_io_elements = n_compute_elements;
	if (MPI_SUCCESS != MPI_Allreduce((const void *)&n_io_elements,
	                                 (void *)&n_io_elements_global,
	                                 1, dtype, MPI_SUM, comm)) {
		return SMIOL_MPI_ERROR;
	}

	/*
	 * Determine the contiguous range of elements to be read/written by
	 * this MPI task
	 */
	ierr = get_io_elements(context->comm_rank, num_io_tasks, io_stride,
	                       n_io_elements_global, &io_start, &io_count);

	/*
	 * Fill in io_elements from io_start through io_start + io_count - 1
	 */
	io_elements = NULL;
	if (io_count > 0) {
		io_elements = (SMIOL_Offset *)malloc(sizeof(SMIOL_Offset)
		                                     * n_io_elements_global);
		if (io_elements == NULL) {
			return SMIOL_MALLOC_FAILURE;
		}
		for (i = 0; i < io_count; i++) {
			io_elements[i] = (SMIOL_Offset)(io_start + i);
		}
	}

	/*
	 * Build the mapping between compute tasks and I/O tasks
	 */
	ierr = build_exchange(context,
	                      n_compute_elements, compute_elements,
	                      io_count, io_elements,
	                      decomp);

	free(io_elements);

	/*
	 * If decomp was successfully created, add io_start and io_count values
	 * to the decomp before returning
	 */
	if (ierr == SMIOL_SUCCESS) {
		(*decomp)->io_start = io_start;
		(*decomp)->io_count = io_count;
	}

	return ierr;
}


/********************************************************************************
 *
 * SMIOL_free_decomp
 *
 * Frees a mapping between compute elements and I/O elements.
 *
 * Free all memory of a SMIOL_decomp and returns SMIOL_SUCCESS. If decomp
 * points to NULL, then do nothing and return SMIOL_SUCCESS. After this routine
 * is called, no other SMIOL routines should use the freed SMIOL_decomp.
 *
 ********************************************************************************/
int SMIOL_free_decomp(struct SMIOL_decomp **decomp)
{
	if ((*decomp) == NULL) {
		return SMIOL_SUCCESS;
	}

	free((*decomp)->comp_list);
	free((*decomp)->io_list);
	free((*decomp));
	*decomp = NULL;

	return SMIOL_SUCCESS;
}


/********************************************************************************
 *
 * build_start_count
 *
 * Constructs start[] and count[] arrays for parallel I/O operations
 *
 * Given a pointer to a SMIOL file that was previously opened, the name of
 * a variable in that file, and a SMIOL decomp, this function returns three
 * items that may be used when reading or writing the variable in parallel:
 *
 * 1) The size of each "element" of the variable, where an element is defined as
 *    a contiguous memory range associated with the slowest-varying, non-record
 *    dimension of the variable; for example, a variable
 *    float foo[nCells][nVertLevels] would have an element size of
 *    sizeof(float) * nVertLevels if nCells were a decomposed dimension.
 *
 * 2) The number of dimensions for the variable, including any unlimited/record
 *    dimension.
 *
 * 3) The start[] and count[] arrays (each with size ndims) to be read or written
 *    by an MPI rank using the I/O decomposition described in decomp.
 *
 * If the decomp argument is NULL, the variable is to be read or written as
 * a non-decomposed variable; typically, only MPI rank 0 will write
 * the non-decomposed variable, and all MPI ranks will read the non-decomposed
 * variable.
 *
 * Depending on the value of the write_or_read argument -- either START_COUNT_READ
 * or START_COUNT_WRITE -- the count[] values will be set so that all ranks will
 * read the variable, or only rank 0 will write the variable.
 *
 ********************************************************************************/
int build_start_count(struct SMIOL_file *file, const char *varname,
                      const struct SMIOL_decomp *decomp,
                      int write_or_read, size_t *element_size, int *ndims,
                      size_t **start, size_t **count)
{
	int i;
	int ierr;
	int vartype;
	char **dimnames;
	SMIOL_Offset *dimsizes;
	int has_unlimited_dim = 0;

/* TO DO - define maximum string size, currently assumed to be 64 chars */

	/*
	 * Figure out type of the variable, as well as its dimensions
	 */
	ierr = SMIOL_inquire_var(file, varname, &vartype, ndims, NULL);
	if (ierr != SMIOL_SUCCESS) {
		return ierr;
	}

	dimnames = malloc(sizeof(char *) * (size_t)(*ndims));
        if (dimnames == NULL) {
		ierr = SMIOL_MALLOC_FAILURE;
		return ierr;
	}

	for (i = 0; i < *ndims; i++) {
		dimnames[i] = malloc(sizeof(char) * (size_t)64);
	        if (dimnames[i] == NULL) {
			int j;

			for (j = 0; j < i; j++) {
				free(dimnames[j]);
			}
			free(dimnames);

			ierr = SMIOL_MALLOC_FAILURE;
			return ierr;
		}
	}

	ierr = SMIOL_inquire_var(file, varname, NULL, NULL, dimnames);
	if (ierr != SMIOL_SUCCESS) {
		for (i = 0; i < *ndims; i++) {
			free(dimnames[i]);
		}
		free(dimnames);
		return ierr;
	}
	
	dimsizes = malloc(sizeof(SMIOL_Offset) * (size_t)(*ndims));
        if (dimsizes == NULL) {
		ierr = SMIOL_MALLOC_FAILURE;
		return ierr;
	}

	/*
	 * It is assumed that only the first dimension can be an unlimited
	 * dimension, so by inquiring about dimensions from last to first, we
	 * can be guaranteed that has_unlimited_dim will be set correctly at
	 * the end of the loop over dimensions
	 */
	for (i = (*ndims-1); i >= 0; i--) {
		ierr = SMIOL_inquire_dim(file, dimnames[i], &dimsizes[i],
		                         &has_unlimited_dim);
		if (ierr != SMIOL_SUCCESS) {
			for (i = 0; i < *ndims; i++) {
				free(dimnames[i]);
			}
			free(dimnames);
			free(dimsizes);

			return ierr;
		}
	}

	for (i = 0; i < *ndims; i++) {
		free(dimnames[i]);
	}
	free(dimnames);

	/*
	 * Set basic size of each element in the field; only necessary if
	 * the field is decomposed and therefore must be transferred prior to
	 * writing
	 */
	*element_size = 1;
	if (decomp) {
		switch (vartype) {
			case SMIOL_REAL32:
				*element_size = sizeof(float);
				break;
			case SMIOL_REAL64:
				*element_size = sizeof(double);
				break;
			case SMIOL_INT32:
				*element_size = sizeof(int);
				break;
			case SMIOL_CHAR:
				*element_size = sizeof(char);
				break;
		}
	}

	*start = malloc(sizeof(size_t) * (size_t)(*ndims));
        if (*start == NULL) {
		free(dimsizes);
		ierr = SMIOL_MALLOC_FAILURE;
		return ierr;
	}

	*count = malloc(sizeof(size_t) * (size_t)(*ndims));
        if (*count == NULL) {
		free(dimsizes);
		free(start);
		ierr = SMIOL_MALLOC_FAILURE;
		return ierr;
	}

	/*
	 * Build start/count description of the part of the variable to be
	 * read or written. Simultaneously, compute the product of all
	 * non-unlimited, non-decomposed dimension sizes, scaled by the basic
	 * element size to get the effective size of each element to be read or
	 * written
	 */
	for (i = 0; i < *ndims; i++) {
		(*start)[i] = (size_t)0;
		(*count)[i] = (size_t)dimsizes[i];

		/*
		 * If variable has an unlimited dimension, set start to current
		 * frame and count to one
		 */
		if (has_unlimited_dim && i == 0) {
			(*start)[i] = (size_t)file->frame;
			(*count)[i] = (size_t)1;
		}

		/*
		 * If variable is decomposed, set the slowest-varying,
		 * non-record dimension start and count based on values from
		 * the decomp structure
		 */
		if (decomp) {
			if ((!has_unlimited_dim && i == 0) ||
			    (has_unlimited_dim && i == 1)) {
				(*start)[i] = decomp->io_start;
				(*count)[i] = decomp->io_count;
			} else {
				*element_size *= (*count)[i];
			}
		} else {
			*element_size *= (*count)[i];
		}

		if (write_or_read == START_COUNT_WRITE) {
			/*
			 * If the variable is not decomposed, only MPI rank 0
			 * will have non-zero count values so that all MPI ranks
			 * do no try to write the same offsets
			 */
			if (!decomp && file->context->comm_rank != 0) {
				(*count)[i] = 0;
			}
		}
	}

	free(dimsizes);

	return SMIOL_SUCCESS;
}
