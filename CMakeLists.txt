# This file is part of enblend/enfuse.
# Licence details can be found in the file COPYING.
#
# Copyright (c) 2009, Kornel Benko <Kornel.Benko@berlin.de>
#                   , Ryan Sleevi <ryan+hugin@sleevi.com>
#                   , Harry van der Wolf <hvdwolf@gmail.com>
# Original file from lyx-project, heavily modified
#
cmake_minimum_required(VERSION 2.4)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
  cmake_policy(SET CMP0005 OLD)
endif(COMMAND cmake_policy)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

project(enblend)
set(CMAKE_BUILD_TYPE Release CACHE STRING
	"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	FORCE)

# where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/
# is checked
# For now, reuse the hugin modules directory
#
# Try to find ...
set(CMAKE_MODULE_PATH) # overwrite cache
find_file(CMAKE_MODULE_PATH CMakeModules
  "${CMAKE_SOURCE_DIR}"
  "${CMAKE_SOURCE_DIR}/../hugin"
  "${CMAKE_SOURCE_DIR}/../hugin-trunk"
  "${CMAKE_SOURCE_DIR}/../../hugin/hugin"
  "${CMAKE_SOURCE_DIR}/../../hugin/trunk"
  "$ENV{HOME}/hugin")
message(STATUS "found hugin modules directory in ${CMAKE_MODULE_PATH}")

include(HuginMacros)
include(FindLibraryForCPU)

set(EXECUTABLE_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/bin)

# Regenerate when calling make
set(CMAKE_SUPPRESS_REGENERATION FALSE)

set(VERSION_INFO "CMake Build")

macro(setstripped _varname _value)
  string(STRIP ${_value} _v)
  if(USE_POSIX_PACKAGING)
    string(TOLOWER ${_v} ${_varname})
  else()
    set(${_varname} ${_v})
  endif()
endmacro(setstripped)

# Try to get the package version from file
set(TOP_SRC_DIR ${CMAKE_SOURCE_DIR})
file(STRINGS "${TOP_SRC_DIR}/VERSION" _version_lines)
foreach(_v_l ${_version_lines})
  #if (_v_l MATCHES "^\([0-9]+\)\.\([0-9]+\)-\([0-9a-z]+\)")
  if (_v_l MATCHES "^\([0-9]+\)\\.\([0-9]+\)-\([0-9a-z]+\)")
    setstripped(ENBLEND_MAJOR_VERSION ${CMAKE_MATCH_1})
    setstripped(ENBLEND_MINOR_VERSION ${CMAKE_MATCH_2})
    message(STATUS "Extracted enblend version = ${ENBLEND_MAJOR_VERSION}.${ENBLEND_MINOR_VERSION}")
    set(PACKAGE_VERSION "${ENBLEND_MAJOR_VERSION}.${ENBLEND_MINOR_VERSION}-hg")
  endif()
endforeach()

## COPIED FROM HUGIN - Needs to run early on to set SOURCE_BASE_DIR for WIN32
IF(WIN32)
  STRING( REGEX REPLACE "(.*)/[^/]+$" "\\1" work "${CMAKE_SOURCE_DIR}" )
  # create the cache entry
  SET( SOURCE_BASE_DIR ${work} CACHE FILEPATH "parent dir of hugin source root" )
  # wxWidgets
  SET(wxWidgets_ROOT_DIR ${SOURCE_BASE_DIR}/wxWidgets-2.8.10)  
ENDIF(WIN32)

# First include so it picks up (our) config.h
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})

SET(wxWidgets_USE_LIBS base core)
FIND_PACKAGE(wxWidgets REQUIRED)

OPTION(ENABLE_OPENMP "OpenMP Support" OFF)
OPTION(ENABLE_IMAGECACHE "Image Cache Support" ON)
OPTION(ENABLE_GPU "GPU Acceleration Support" OFF)
IF(NOT CMAKE_CL_64)
  OPTION(ENABLE_SSE2 "SSE2 Support" OFF)
ELSE(NOT CMAKE_CL_64)
  SET(ENABLE_SSE2 OFF)
ENDIF(NOT CMAKE_CL_64)

IF(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  add_definitions(${OpenMP_CXX_FLAGS})
ENDIF(ENABLE_OPENMP)

IF(ENABLE_IMAGECACHE AND ENABLE_OPENMP)
  MESSAGE(FATAL_ERROR "The image cache does not currently work with OpenMP. You must disable either the image cache or OpenMP.")
ENDIF(ENABLE_IMAGECACHE AND ENABLE_OPENMP)


IF(ENABLE_IMAGECACHE)
  ADD_DEFINITIONS("-DENBLEND_CACHE_IMAGES")
ENDIF(ENABLE_IMAGECACHE)

IF(ENABLE_GPU)
  find_package(GLUT REQUIRED)
  find_package(GLEW REQUIRED)
  
  INCLUDE_DIRECTORIES(${GLUT_INCLUDE_DIR} ${GLEW_INCLUDE_DIR})
  set(common_libs ${GLUT_LIBRARIES} ${GLEW_LIBRARIES})
  ADD_DEFINITIONS(-DHAVE_LIBGLEW)
ENDIF(ENABLE_GPU)

IF(ENABLE_SSE2)
  IF(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -msse2 -mtune=pentium4")
  ELSEIF(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /arch:SSE2")
  ENDIF(CMAKE_COMPILER_IS_GNUCXX)
ENDIF(ENABLE_SSE2)

IF(MSVC)
  ADD_DEFINITIONS("-D_CRT_SECURE_NO_DEPRECATE")
ENDIF(MSVC)

find_library(LIBZ z)

FIND_PACKAGE(TIFF REQUIRED)
FIND_PACKAGE(JPEG REQUIRED)
FIND_PACKAGE(PNG REQUIRED)
FIND_PACKAGE(OpenEXR REQUIRED)
FIND_PACKAGE(LCMS REQUIRED)
FIND_PACKAGE(LibXMI REQUIRED)
IF(MSVC)
  FIND_PACKAGE(STLport REQUIRED)
ENDIF(MSVC)
FIND_PACKAGE(Threads)
SET(HasTIFF ${TIFF_FOUND})
SET(HasJPEG ${JPEG_FOUND})
SET(HasPNG ${PNG_FOUND})
SET(HasEXR ${OPENEXR_FOUND})

##
## Boost
##
IF(WIN32)
  SET(Boost_USE_STATIC_LIBS ON)
  # If BOOST_ROOT was defined in the environment, use it. 
  if (NOT BOOST_ROOT AND NOT $ENV{BOOST_ROOT} STREQUAL "") 
    set(BOOST_ROOT $ENV{BOOST_ROOT}) 
  else(NOT BOOST_ROOT AND NOT $ENV{BOOST_ROOT} STREQUAL "") 
    set(BOOST_ROOT ${SOURCE_BASE_DIR}/boost_1_39_0) 
   endif(NOT BOOST_ROOT AND NOT $ENV{BOOST_ROOT} STREQUAL "") 
   SET(BOOST_DIR_SEARCH_USER ${BOOST_ROOT})
ENDIF(WIN32)
FIND_PACKAGE(Boost REQUIRED thread date_time)
IF (NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR "boost library not found. If it is installed with a\nunrecognized suffix, specify it with -DBOOST_LIBRARYDIR=<librarydir> and -DBOOST_INCLUDEDIR=<includedir>")
ENDIF (NOT Boost_FOUND)


INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ${STLport_INCLUDE_DIR})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS} ${STLport_LIBRARIES_DIR})

MESSAGE(STATUS "Boost inc dirs = ${Boost_INCLUDE_DIRS}")

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR} ${TIFF_INCLUDE_DIR} ${JPEG_INCLUDE_DIR} ${PNG_INCLUDE_DIR} ${LCMS_INCLUDE_DIR} ${LIBXMI_INCLUDE_DIR} ${ZLIB_INCLUDE_DIR} ${OPENEXR_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${iconv_INCLUDE_DIR})
set(common_libs ${common_libs} ${Boost_LIBRARIES} ${JPEG_LIBRARIES} ${TIFF_LIBRARIES} ${PNG_LIBRARIES} ${LCMS_LIBRARIES} ${LIBXMI_LIBRARIES})

# here I don't know, how to get the revision number from mercurial
# so the try is "hg tip"
FIND_PROGRAM(_hg hg)
if (NOT ${_hg} MATCHES "-NOTFOUND")
  EXECUTE_PROCESS(COMMAND ${_hg} tip WORKING_DIRECTORY "${TOP_SRC_DIR}" OUTPUT_VARIABLE _release OUTPUT_STRIP_TRAILING_WHITESPACE)
  foreach(_v_l ${_release})
    if (_v_l MATCHES "^changeset: *[^0-9]*\([0-9]+\):\([a-z0-9]+\)")
      set(CPACK_PACKAGE_VERSION_PATCH ${CMAKE_MATCH_1})
      set(PACKAGE_REVISION ${CMAKE_MATCH_2})
    endif()
  endforeach()
else()
  set(CPACK_PACKAGE_VERSION_PATCH 468)
endif()

# Try to create config.h
set(ENBLEND_VERSION_ONLY "${ENBLEND_MAJOR_VERSION}.${ENBLEND_MINOR_VERSION}-${PACKAGE_REVISION}")
set(PACKAGE_STRING "enblend-enfuse ${ENBLEND_VERSION_ONLY}")
set(ENBLEND_VERSION "${PACKAGE_STRING}")
include(ConfigureChecks.cmake)
configure_file(config.h.cmake ${CMAKE_BINARY_DIR}/config.h)

set(library_type STATIC)
add_definitions(-DVIGRA_STATIC_LIB)

add_definitions(-DHAVE_CONFIG_H)
include_directories(${CMAKE_BINARY_DIR} ${TOP_SRC_DIR}/include)

IF (APPLE)
   #This namespace workaround is not needed for CMake >= 2.5.x
   ## Add needed Compiler and Linker flags for OSX
      IF (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
         SET (CMAKE_C_FLAGS:STRING "\"$CMAKE_C_FLAGS\" -flat_namespace -undefined suppress" )
         SET (CMAKE_CXX_FLAGS:STRING "\"$CMAKE_CXX_FLAGS\" -flat_namespace -undefined suppress" )
         SET (CMAKE_LDFLAGS:STRING "\"$CMAKE_LD_FLAGS\"" )
         SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  -flat_namespace -undefined suppress" )
         SET (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS}  -flat_namespace -undefined suppress" )
         SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}  -flat_namespace -undefined suppress" )
      ENDIF (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
ENDIF (APPLE)

add_subdirectory(src)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bin/enblend ${CMAKE_CURRENT_BINARY_DIR}/bin/enfuse DESTINATION bin CONFIGURATIONS Release)

SET(CPACK_PACKAGE_VERSION_MAJOR ${ENBLEND_MAJOR_VERSION})
SET(CPACK_PACKAGE_VERSION_MINOR ${ENBLEND_MINOR_VERSION})

# needed by rpm
SET(CPACK_SET_DESTDIR "ON")

# The next ones are required by debian
SET(CPACK_PACKAGE_CONTACT "hugin-ptx@googlegroups.com")

include(CPack)

