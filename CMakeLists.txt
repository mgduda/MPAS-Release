cmake_minimum_required(VERSION 3.0)
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
    add_compile_definitions(CORE_ATMOSPHERE)
elseif(MPAS_CORE STREQUAL "ocean")
    set(MPAS_EXE_NAME "ocean_model")
    add_compile_definitions(CORE_OCEAN)
else()
    message(FATAL_ERROR "Please define MPAS_CORE as one of: atmosphere ocean")
endif()


#
# Set compiler and compiler options
#
if(MPAS_COMPILER STREQUAL gnu)
    set(CMAKE_C_COMPILER mpicc)
    set(CMAKE_Fortran_COMPILER mpifort)
elseif(MPAS_COMPILER STREQUAL intel)
    set(CMAKE_C_COMPILER icc)
    set(CMAKE_Fortran_COMPILER ifort)
else()
    message(FATAL_ERROR "Please define MPAS_COMPILER as one of: gnu intel")
endif()


#
# Add required libraries: netCDF, parallel-netCDF, PIO
#
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH})

find_package(NetCDF REQUIRED)
find_package(PNETCDF REQUIRED)

#
# Define main target and dependencies
#
add_executable(${MPAS_EXE_NAME})

include_directories(${PNETCDF_INCLUDE_DIRS} ${NetCDF_INCLUDE_DIRS})
target_link_libraries(${MPAS_EXE_NAME} ${NetCDF_LIBRARIES} ${PNETCDF_LIBRARIES} ${CMAKE_DL_LIBS})

set_property(TARGET ${MPAS_EXE_NAME} PROPERTY Fortran_FORMAT FREE)

add_subdirectory(src)
