# - Find TCMALLOC
# Find the native TCMALLOC includes and library
#
#  TCMALLOC_INCLUDE_DIR - where to find tcmalloc.h etc.
#  TCMALLOC_LIBRARIES   - List of libraries when using tcmalloc.
#  TCMALLOC_FOUND       - True if tcmalloc found.

FIND_PATH(TCMALLOC_INCLUDE_DIR tcmalloc.h
  PATHS
    /usr/local/include
    /usr/include
    ${CMAKE_SYSTEM_INCLUDE_PATH}/gperftools
    ${CMAKE_INCLUDE_PATH}/gperftools
    ${SOURCE_BASE_DIR}/gperftools-2.0/src/windows/gperftools
)

find_library(TCMALLOC_LIBRARIES 
  NAMES tcmalloc tcmalloc_minimal libtcmalloc_minimal
  PATHS 
    ${SYSTEM_LIB_DIRS}
    ${SOURCE_BASE_DIR}/gperftools-2.0/x64/Release    
    ${SOURCE_BASE_DIR}/gperftools-2.0/Release    
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TCMALLOC DEFAULT_MSG 
                                  TCMALLOC_INCLUDE_DIR TCMALLOC_LIBRARIES)

MARK_AS_ADVANCED(
  TCMALLOC_LIBRARIES
  TCMALLOC_INCLUDE_DIR
  )
