# - Try to find netCDF
# Once done this will define
#  NetCDF_FOUND - System has netCDF
#  NetCDF_INCLUDE_DIRS - The netCDF include directories
#  NetCDF_LIBRARIES - The libraries needed to use netCDF

find_path(NETCDF_INCLUDE_DIR netcdf.mod NETCDF.mod HINTS ENV NETCDF ENV NETCDF_PATH PATH_SUFFIXES include)
find_library(NETCDF_Fortran_LIBRARY libnetcdff.a HINTS ENV NETCDF ENV NETCDF_PATH PATH_SUFFIXES lib)
find_library(NETCDF_C_LIBRARY libnetcdf.a HINTS ENV NETCDF ENV NETCDF_PATH PATH_SUFFIXES lib)

set(NetCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR} )
set(NetCDF_LIBRARIES ${NETCDF_Fortran_LIBRARY} ${NETCDF_C_LIBRARY})

find_library(HDF5HL_LIBRARY libhdf5_hl.a HINTS ENV NETCDF ENV NETCDF_PATH PATH_SUFFIXES lib)
if(HDF5HL_LIBRARY)
    set(NetCDF_LIBRARIES ${NetCDF_LIBRARIES} ${HDF5HL_LIBRARY})
else()
    message("HDF5 HL library not found...")
endif()

find_library(HDF5_LIBRARY libhdf5.a HINTS ENV NETCDF ENV NETCDF_PATH PATH_SUFFIXES lib)
if(HDF5_LIBRARY)
    set(NetCDF_LIBRARIES ${NetCDF_LIBRARIES} ${HDF5_LIBRARY})
else()
    message("HDF5 library not found...")
endif()

find_library(ZLIB_LIBRARY libz.a HINTS ENV NETCDF ENV NETCDF_PATH PATH_SUFFIXES lib)
if(ZLIB_LIBRARY)
    set(NetCDF_LIBRARIES ${NetCDF_LIBRARIES} ${ZLIB_LIBRARY})
else()
    message("ZLIB library not found...")
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(NetCDF DEFAULT_MSG NetCDF_LIBRARIES NetCDF_INCLUDE_DIRS)
#find_package_handle_standard_args(NetCDF "Try setting the NETCDF environment variable" NetCDF_LIBRARIES NetCDF_INCLUDE_DIRS)
