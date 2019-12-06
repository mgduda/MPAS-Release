# - Try to find PIO
# Once done this will define
#  PIO_FOUND - System has netCDF
#  PIO_INCLUDE_DIRS - The netCDF include directories
#  PIO_LIBRARIES - The libraries needed to use netCDF
#  PIO_VERSION - A global property indicating the PIO API: either 1 or 2

find_path(PIO_INCLUDE_DIR pio.mod PIO.mod HINTS ENV PIO ENV PIO_PATH PATH_SUFFIXES include)
find_library(PIO1_C_LIBRARY libpio.a HINTS ENV PIO ENV PIO_PATH PATH_SUFFIXES lib)
find_library(PIO2_C_LIBRARY libpioc.a HINTS ENV PIO ENV PIO_PATH PATH_SUFFIXES lib)
find_library(PIO_F_LIBRARY libpiof.a HINTS ENV PIO ENV PIO_PATH PATH_SUFFIXES lib)

if(PIO1_C_LIBRARY)
   set(PIO_C_LIBRARY ${PIO1_C_LIBRARY})
   set_property(GLOBAL PROPERTY PIO_VERSION 1)
elseif(PIO2_C_LIBRARY)
   set(PIO_C_LIBRARY ${PIO2_C_LIBRARY})
   set_property(GLOBAL PROPERTY PIO_VERSION 2)
endif()

set(PIO_INCLUDE_DIRS ${PIO_INCLUDE_DIR})
set(PIO_LIBRARIES ${PIO_F_LIBRARY} ${PIO_C_LIBRARY})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PIO_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(PIO DEFAULT_MSG PIO_LIBRARIES PIO_INCLUDE_DIRS)
