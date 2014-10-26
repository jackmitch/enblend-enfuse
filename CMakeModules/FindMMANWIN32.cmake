# - Find MMANWIN32
# Find the wrapper around sys/mman.h for Windows
#
#  MMAN_WIN32_INCLUDE_DIR - where to find mman.h etc.
#  MMAN_WIN32_LIBRARIES   - List of libraries when using tcmalloc.
#  MMAN_WIN32_FOUND       - True if tcmalloc found.

FIND_PATH(MMAN_WIN32_INCLUDE_DIR mman.h
  PATHS
    ${SOURCE_BASE_DIR}/mman-win32
)

find_library(MMAN_WIN32_LIBRARIES
  NAMES mman.lib
  PATHS
    ${SYSTEM_LIB_DIRS}
    ${SOURCE_BASE_DIR}/mman-win32/x64/Release
    ${SOURCE_BASE_DIR}/mman-win32/Release
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MMAN_WIN32 DEFAULT_MSG
                                  MMAN_WIN32_INCLUDE_DIR MMAN-MMAN_WIN32_LIBRARIES)

MARK_AS_ADVANCED(
  MMAN_WIN32_LIBRARIES
  MMAN-WIN32_INCLUDE_DIR
)
