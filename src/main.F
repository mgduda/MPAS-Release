program hello

    use netcdf, only : nf90_create, nf90_def_dim, nf90_enddef, nf90_close, nf90_clobber, nf90_unlimited
    use pnetcdf, only : nf90mpi_create, nf90mpi_def_dim, nf90mpi_enddef, nf90mpi_close
    use mpi
    use pio
    use pio_types
    use ESMF
    use ESMF_ClockMod

#ifdef CORE_ATMOSPHERE
    use atm_core
#endif


    implicit none

    interface
       subroutine parse_xml_streams() bind(c)
       end subroutine parse_xml_streams
    end interface

#ifdef SINGLE_PRECISION
    integer, parameter :: RKIND  = selected_real_kind(6)
#else
    integer, parameter :: RKIND  = selected_real_kind(12)
#endif

    integer :: ierr, ncid, dimid
    integer (kind=MPI_OFFSET_KIND) :: timeL = 1
    type (iosystem_desc_t) :: io_system
    type (ESMF_Clock) :: clock
    type (ESMF_Time) :: time

    real(kind=RKIND), parameter :: pi = 4.0_RKIND * atan(1.0_RKIND)
    real :: sqrt_two

#ifdef CORE_ATMOSPHERE
    write(0,*) 'This is MPAS-Atmosphere'
    call atm_setup_core()
    call atm_setup_domain()
#endif
#ifdef CORE_OCEAN
    write(0,*) 'This is MPAS-Ocean'
#endif

    write(0,*) 'pi = ', pi
    sqrt_two = 2.0
    sqrt_two = sqrt(sqrt_two)
    write(0,*) 'sqrt(2.0) = ', sqrt_two

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

    write(0,*) 'Here is a long source line to test whether the appropriate compiler flags have been added to permit source lines longer than 132 characters'

    stop

end program hello
