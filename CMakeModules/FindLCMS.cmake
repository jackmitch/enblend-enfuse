
IF(NOT WIN32)
	FIND_LIBRARY(LCMS_LIBRARIES lcms)
ELSE(NOT WIN32)
	FIND_PATH(LCMS_ROOT_DIR
	  NAMES include/lcms.h
	  PATHS /usr/local
	    /usr
	    ${SOURCE_BASE_DIR}
	  PATH_SUFFIXES
	    lcms-1.17
		lcms-1.18
	)
	
	FIND_PATH(LCMS_INCLUDE_DIR 
	  NAMES lcms.h
	  PATHS
	    /usr/local/include
	    /usr/include
	    ${LCMS_ROOT_DIR}/include
	)

	include(FindLibraryWithDebug)
	find_library_with_debug(LCMS_LIBRARIES
	  WIN32_DEBUG_POSTFIX d	
	  NAMES lcms
	  PATHS ${LCMS_ROOT_DIR}/Lib/MS
	)
	
	MARK_AS_ADVANCED(
	  LCMS_ROOT_DIR
	  LCMS_LIBRARIES
	  LCMS_INCLUDE_DIR
	  )		
ENDIF(NOT WIN32) 