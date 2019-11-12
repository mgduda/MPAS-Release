program hello

    use netcdf
    use mpi

    integer :: ierr, ncid
#ifdef CORE_ATMOSPHERE
    write(0,*) 'This is MPAS-Atmosphere'
#endif
#ifdef CORE_OCEAN
    write(0,*) 'This is MPAS-Ocean'
#endif

    call MPI_Init(ierr)
    ierr = nf90_open('blah.nc', NF90_CLOBBER, ncid)

    ierr = nf90_close(ncid)
    call MPI_Finalize(ierr)

    stop

end program hello
