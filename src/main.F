program hello

    use netcdf, only : nf90_create, nf90_def_dim, nf90_enddef, nf90_close, nf90_clobber, nf90_unlimited
    use pnetcdf, only : nf90mpi_create, nf90mpi_def_dim, nf90mpi_enddef, nf90mpi_close
    use mpi
    use pio
    use pio_types

    implicit none

    interface
       subroutine parse_xml_streams() bind(c)
       end subroutine parse_xml_streams
    end interface

    integer :: ierr, ncid, dimid
    integer (kind=MPI_OFFSET_KIND) :: timeL = 1
    type (iosystem_desc_t) :: io_system

#ifdef CORE_ATMOSPHERE
    write(0,*) 'This is MPAS-Atmosphere'
#endif
#ifdef CORE_OCEAN
    write(0,*) 'This is MPAS-Ocean'
#endif

    call MPI_Init(ierr)

    ! NetCDF
    ierr = nf90_create('blah.nc', NF90_CLOBBER, ncid)
    ierr = nf90_def_dim(ncid, 'Time', NF90_UNLIMITED, dimid)
    ierr = nf90_enddef(ncid)
    ierr = nf90_close(ncid)

    ! PNetCDF
    ierr = nf90mpi_create(MPI_COMM_WORLD, 'pnetcdfBlash.nc', NF90_CLOBBER, MPI_INFO_NULL, ncid)
    ierr = nf90mpi_def_dim(ncid, 'Time', timeL, dimid)
    ierr = nf90mpi_enddef(ncid)
    ierr = nf90mpi_close(ncid)

    ! PIO
    call PIO_init(0, MPI_COMM_WORLD, 1, 0, 1, PIO_rearr_box, io_system)

    call MPI_Finalize(ierr)

    call parse_xml_streams()

    stop

end program hello
