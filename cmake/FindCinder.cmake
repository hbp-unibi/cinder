# - Try to find cinder
# Once done this will define
#  CINDER_FOUND - System has cinder
#  CINDER_INCLUDE_DIRS - The cinder include directories
#  CINDER_LIBRARIES - The libraries needed to use cinder

find_path(CINDER_INCLUDE_DIR cinder/config.h PATH_SUFFIXES cinder)
find_library(CINDER_LIBRARY NAMES cinder libcinder)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cinder DEFAULT_MSG
                                  CINDER_LIBRARY CINDER_INCLUDE_DIR)

mark_as_advanced(CINDER_INCLUDE_DIR CINDER_LIBRARY)

set(CINDER_LIBRARIES ${CINDER_LIBRARY})
set(CINDER_INCLUDE_DIRS ${CINDER_INCLUDE_DIR})
