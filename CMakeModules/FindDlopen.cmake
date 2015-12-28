# - Find DLOPEN
# Find the native DLOPEN includes and library
#
#  DLOPEN_INCLUDE_DIR - where to find dlfcn.h etc.
#  DLOPEN_LIBRARIES   - List of libraries when using dlopen.
#  DLOPEN_FOUND       - True if dlopen found.

FIND_PATH(DLOPEN_INCLUDE_DIR dlfcn.h
  PATHS
    /usr/local/include
    /usr/include
)

find_library(DLOPEN_LIBRARIES 
  NAMES dlopen dl
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DLOPEN DEFAULT_MSG 
                                  DLOPEN_INCLUDE_DIR DLOPEN_LIBRARIES)

MARK_AS_ADVANCED(
  DLOPEN_LIBRARIES
  DLOPEN_INCLUDE_DIR
  )
