# - Find Vistle
# Find the native lvistle_sensei library
#
# set VISTLE_DIR to the build directory
# For example on my workstation, I use
# -DVISTLE_DIR=~/vistle/build-linux64-ompi 
# lib
if(VISTLE_INCLUDE_DIR AND VISTLE_LIBRARY)
  set(VISTLE_FIND_QUIETLY TRUE)
endif()

if(VISTLE_DIR)
  find_library(VISTLE_LIBRARY NAMES libvistle_sensei
    PATHS "${VISTLE_DIR}/lib"
    NO_DEFAULT_PATH)
  set(VISTLE_INCLUDE_DIR "${VISTLE_DIR}/../vistle/insitu/sensei")
    
endif()

find_library(VISTLE_LIBRARY NAMES libvistle_sensei
  PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH)

find_library(VISTLE_LIBRARY NAMES libvistle_sensei)

if(VISTLE_LIBRARY)
  get_filename_component(VISTLE_LIB_DIR ${VISTLE_LIBRARY} DIRECTORY)
  set(VISTLE_INCLUDE_DIR "${VISTLE_LIB_DIR}/../../vistle/insitu/sensei")
  mark_as_advanced(VISTLE_LIBRARY)
  mark_as_advanced(VISTLE_INCLUDE_DIR)
else()
    message("VISTLE_LIBRARY is empty")
endif()

set(err_msg "Failed to locate vistle_sensei in VISTLE_DIR=\"${VISTLE_DIR}\"")

# validate
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VISTLE
  ${err_msg} VISTLE_LIBRARY VISTLE_INCLUDE_DIR)
