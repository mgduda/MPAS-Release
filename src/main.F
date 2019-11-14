program hello

    use netcdf
    use mpi

    interface
       subroutine parse_xml_streams() bind(c)
       end subroutine parse_xml_streams
    end interface

    integer :: ierr, ncid
#ifdef CORE_ATMOSPHERE
    write(0,*) 'This is MPAS-Atmosphere'
#endif
#ifdef CORE_OCEAN
    write(0,*) 'This is MPAS-Ocean'
#endif

    call MPI_Init(ierr)
    ! NetCDF
    ierr = nf90_create('blah.nc', NCF90_CLOBBER, ncid)
    ierr = nf90_def_dim(ncid, "Time", 1, dimid)
    ierr = nf90_enddef(ncid)
    ierr = nf90_close(ncid)
    call MPI_Finalize(ierr)

    call parse_xml_streams()

    stop

end program hello
