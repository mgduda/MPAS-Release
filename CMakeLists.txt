cmake_minimum_required(VERSION 3.11)
cmake_policy(SET CMP0076 NEW)
project(reg C Fortran)



#
# Define build options
#
set(MPAS_CORE "None" CACHE STRING "The MPAS core to build. Possible options: atmosphere ocean")
set(MPAS_COMPILER "None" CACHE STRING "The MPAS core to build. Possible options: gnu intel")
set(MPAS_PRECISION "double" CACHE STRING "The default precision for real values. Possible options: single double")


#
# Set core to build and associated preprocessing definitions
#
if(MPAS_CORE STREQUAL "atmosphere")
    set(MPAS_EXE_NAME "atmosphere_model")
    set(MPAS_CORE_DEF "CORE_ATMOSPHERE")
elseif(MPAS_CORE STREQUAL "ocean")
    set(MPAS_EXE_NAME "ocean_model")
    set(MPAS_CORE_DEF "CORE_OCEAN")
else()
    message(FATAL_ERROR "Please define MPAS_CORE as one of: atmosphere ocean")
endif()

add_executable(${MPAS_EXE_NAME})
target_compile_definitions(${MPAS_EXE_NAME} PUBLIC ${MPAS_CORE_DEF})


#
# Set compilation macros that are invariant
#
target_compile_definitions(${MPAS_EXE_NAME} PUBLIC _MPI)


#
# Set compiler and compiler options
#
if(MPAS_COMPILER STREQUAL gnu)
    set(CMAKE_C_COMPILER mpicc)
    set(CMAKE_Fortran_COMPILER mpif90)
    set(CMAKE_Fortran_FLAGS "-ffree-line-length-none")
    set(PROMOTION_FLAGS "-fdefault-real-8 -fdefault-double-8")
elseif(MPAS_COMPILER STREQUAL intel)
    set(CMAKE_C_COMPILER mpicc)
    set(CMAKE_Fortran_COMPILER mpif90)
    set(PROMOTION_FLAGS "-real-size 64")
else()
    message(FATAL_ERROR "Please define MPAS_COMPILER as one of: gnu intel")
endif()


#
# Set precision for real values
#
if(MPAS_PRECISION STREQUAL "single")
    target_compile_definitions(${MPAS_EXE_NAME} PUBLIC SINGLE_PRECISION)
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${PROMOTION_FLAGS}")
endif()


#
# Add required external libraries: netCDF, parallel-netCDF, PIO
#
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH})

find_package(NetCDF REQUIRED)
find_package(PNETCDF REQUIRED)
find_package(PIO REQUIRED)


#
# Set compilation macros to select appropriate PIO API
#
get_property(PIO_VERSION GLOBAL PROPERTY PIO_VERSION)
if(PIO_VERSION STREQUAL "2")
    target_compile_definitions(${MPAS_EXE_NAME} PUBLIC USE_PIO2)
endif()


#
# Define main target dependencies
#
target_include_directories(${MPAS_EXE_NAME}
                           PUBLIC
                              ${PNETCDF_INCLUDE_DIRS}
                              ${PIO_INCLUDE_DIR}
                              ${NETCDF_INCLUDE_DIRS}
                          )

target_link_libraries(${MPAS_EXE_NAME}
                      PUBLIC
                         ${PIO_LIBRARIES}
                         ${NetCDF_LIBRARIES}
                         ${PNETCDF_LIBRARIES}
                         ${CMAKE_DL_LIBS}
                     )

set_property(TARGET ${MPAS_EXE_NAME} PROPERTY Fortran_FORMAT FREE)


#
# Process everything in the src/ subdirectory
#
add_subdirectory(src)


#
# Add include paths and libraries from packages that are internally built by MPAS
#
get_property(MPAS_INTERNAL_INCLUDE_DIRS GLOBAL PROPERTY MPAS_INTERNAL_INCLUDE_DIRS)
get_property(MPAS_INTERNAL_LIBRARIES GLOBAL PROPERTY MPAS_INTERNAL_LIBRARIES)
target_include_directories(${MPAS_EXE_NAME} PRIVATE ${MPAS_INTERNAL_INCLUDE_DIRS})
target_link_libraries(${MPAS_EXE_NAME} PRIVATE ${MPAS_INTERNAL_LIBRARIES})
