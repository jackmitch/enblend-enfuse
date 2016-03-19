IF(NOT WIN32)
  FIND_PATH(LCMS2_INCLUDE_DIR
    NAMES lcms2.h
    PATHS 
      /usr/local/include
      /usr/include
  )
  FIND_LIBRARY(LCMS2_LIBRARIES lcms2 HINTS /usr/local/lib /usr/lib/x86_64-linux-gnu /usr/lib32)
ELSE(NOT WIN32)
    FIND_PATH(LCMS2_ROOT_DIR
      NAMES include/lcms2.h
      PATHS /usr/local
        /usr
        ${SOURCE_BASE_DIR}
      PATH_SUFFIXES
        lcms2-2.7
        lcms2-2.6
        lcms2-2.5
    )
    
    FIND_PATH(LCMS2_INCLUDE_DIR 
      NAMES lcms2.h
      PATHS ${LCMS2_ROOT_DIR}/include
    )

    include(FindLibraryWithDebug)
    find_library_with_debug(LCMS2_LIBRARIES
      WIN32_DEBUG_POSTFIX d    
      NAMES lcms2 lcms2_static
      PATHS ${LCMS2_ROOT_DIR}/Lib/MS
    )
    
    MARK_AS_ADVANCED(LCMS2_ROOT_DIR)
ENDIF(NOT WIN32) 

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LCMS2 DEFAULT_MSG 
                                  LCMS2_INCLUDE_DIR LCMS2_LIBRARIES)
MARK_AS_ADVANCED(LCMS2_LIBRARIES LCMS2_INCLUDE_DIR)
