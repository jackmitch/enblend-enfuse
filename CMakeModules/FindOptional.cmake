IF (OPTONAL_INCLUDE_DIR)
  # Already in cache, be silent
  SET(OPTOINAL_FIND_QUIETLY TRUE)
ENDIF()

FIND_PATH(OPTIONAL_INCLUDE_DIR optional.hpp
  /usr/local/include
  /usr/include
  ${SOURCE_BASE_DIR}/Optional-master
)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Optional DEFAULT_MSG 
                                  OPTIONAL_INCLUDE_DIR)

MARK_AS_ADVANCED(
  OPTIONAL_INCLUDE_DIR
)
