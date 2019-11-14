add_subdirectory(ezxml)
add_subdirectory(esmf_time_f90)

target_include_directories(${MPAS_EXE_NAME} PUBLIC esmf_time_f90)
