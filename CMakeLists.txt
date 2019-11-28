cmake_minimum_required(VERSION 3.5)
cmake_policy(SET CMP0076 NEW)
project(reg C Fortran)



#
# Define build options
#
set(MPAS_CORE "None" CACHE STRING "The MPAS core to build. Possible options: atmosphere ocean")
set(MPAS_COMPILER "None" CACHE STRING "The MPAS core to build. Possible options: gnu intel")


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

################################################################################
##### Temporary
################################################################################
target_compile_definitions(${MPAS_EXE_NAME} PUBLIC USE_PIO2)
target_compile_definitions(${MPAS_EXE_NAME} PUBLIC _MPI)
target_compile_definitions(${MPAS_EXE_NAME} PUBLIC DO_PHYSICS)
target_compile_definitions(${MPAS_EXE_NAME} PUBLIC mpas)
target_compile_definitions(${MPAS_EXE_NAME} PUBLIC SINGLE_PRECISION)
target_compile_definitions(${MPAS_EXE_NAME} PUBLIC MPAS_NAMELIST_SUFFIX=atmosphere)


#
# Set compiler and compiler options
#
if(MPAS_COMPILER STREQUAL gnu)
    set(CMAKE_C_COMPILER mpicc)
    set(CMAKE_Fortran_COMPILER mpifort)
    set(CMAKE_Fortran_FLAGS "-ffree-line-length-none")
elseif(MPAS_COMPILER STREQUAL intel)
    set(CMAKE_C_COMPILER icc)
    set(CMAKE_Fortran_COMPILER ifort)
else()
    message(FATAL_ERROR "Please define MPAS_COMPILER as one of: gnu intel")
endif()


#
# Add required external libraries: netCDF, parallel-netCDF, PIO
#
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH})

find_package(NetCDF REQUIRED)
find_package(PNETCDF REQUIRED)
find_package(PIO REQUIRED)


#
# Define main target dependencies
#
target_include_directories(${MPAS_EXE_NAME}
                           PUBLIC
                              ${PNETCDF_INCLUDE_DIRS}
                              ${NETCDF_INCLUDE_DIRS}
                              ${PIO_INCLUDE_DIR})

target_link_libraries(${MPAS_EXE_NAME}
                      PUBLIC
                         ${NetCDF_LIBRARIES}
                         ${PNETCDF_LIBRARIES}
                         ${PIO_LIBRARIES}
                         ${CMAKE_DL_LIBS})

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
