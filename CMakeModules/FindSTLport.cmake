# - Find the STLport includes and libraries.
# The following variables are set if STLport is found.  If STLport is not
# found, STLport_FOUND is set to false.
#  STLport_FOUND                  - True when the STLport include directory is found.
#  STLport_INCLUDE_DIR            - the path to where the STLport include files are.
#  STLport_LIBRARIES_DIR          - the path to where the STLport libraries are
#
# TODO: Support other platforms than Win32
IF(MSVC)
  SET(SUFFIX_FOR_PATH
    STLport-5.2.1
	STLport-trunk
  )
  
  SET(STLPORT_DIR_SEARCH
    ${SOURCE_BASE_DIR}
  )
  
  FIND_PATH(STLport_ROOT_DIR
    NAMES stlport/assert.h
	PATHS
	  ${STLPORT_DIR_SEARCH}
	PATH_SUFFIXES
	  ${SUFFIX_FOR_PATH}
  )
  
  FIND_PATH(STLport_INCLUDE_DIR
    NAMES assert.h
	PATHS ${STLport_ROOT_DIR}/stlport
  )
  
  # TODO: This isn't right... need to adopt a boost-style magic naming depending on what they're linking
  SET(STLport_LIBRARIES_DIR ${STLport_ROOT_DIR}/lib)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(STLport DEFAULT_MSG STLport_INCLUDE_DIR STLport_LIBRARIES_DIR)
  
  MARK_AS_ADVANCED(
    STLport_INCLUDE_DIR
  )		  
ENDIF(MSVC)