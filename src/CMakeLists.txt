# This file is part of enblend/enfuse.
# Licence details can be found in the file COPYING.
#
# Copyright (c) 2009-2013, Kornel Benko <Kornel.Benko@berlin.de>
#                   , Ryan Sleevi <ryan+hugin@sleevi.com>
#                   , Harry van der Wolf <hvdwolf@gmail.com>
#

# Get the version-string from ${TOP_SRC_DIR}/VERSION

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/signature.h"
  COMMAND ${PERL_EXECUTABLE} -I"${TOP_SRC_DIR}/src" "${TOP_SRC_DIR}/src/gen_sig"
  "--extra=${ENBLEND_MAJOR_VERSION}.${ENBLEND_MINOR_VERSION}-${ENBLEND_CHANGESET}"
  ">" "${CMAKE_CURRENT_BINARY_DIR}/signature.h"
  DEPENDS "${TOP_SRC_DIR}/src/gen_sig" "${TOP_SRC_DIR}/VERSION"
  )
set_source_files_properties("${CMAKE_CURRENT_BINARY_DIR}/signature.h" GENERATED)

# prepare all OpenCL sources
SET(OPENCL_SOURCES)
FILE(GLOB opencl_sources ${TOP_SRC_DIR}/src/*.cl)
FOREACH(_cl_source ${opencl_sources})
  GET_FILENAME_COMPONENT(_cl_filename "${_cl_source}" NAME)
  
  STRING(REGEX REPLACE ".cl$" ".icl" _icl_source "${_cl_filename}")
  SET(_icl_source "${CMAKE_CURRENT_BINARY_DIR}/kernels/${_icl_source}")
  STRING(REGEX REPLACE ".cl$" "_source_code" _cl_label "${_cl_filename}")

  # if you want separate OpenCL source files install them now  
  IF(PREFER_SEPARATE_OPENCL_SOURCE)
    INSTALL(FILES ${_cl_source} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/enblend/kernels)
  ELSE()
    ADD_CUSTOM_COMMAND(
      OUTPUT "${_icl_source}"
      COMMAND ${PERL_EXECUTABLE} "${TOP_SRC_DIR}/src/embrace"
        "--format=c++" "--label=${_cl_label}" "${_cl_source}" ">" "${_icl_source}"
      DEPENDS "${_cl_source}"
      COMMENT "Generating ${_icl_source}"
    )
    SET_SOURCE_FILES_PROPERTIES("${_icl_source}" GENERATED)
    LIST(APPEND OPENCL_SOURCES "${_icl_source}")
  ENDIF()
ENDFOREACH()

IF(NOT ${PREFER_SEPARATE_OPENCL_SOURCE})
  ADD_CUSTOM_TARGET(cl_sources ALL DEPENDS "${OPENCL_SOURCES}")
  INCLUDE_DIRECTORIES("${CMAKE_CURRENT_BINARY_DIR}/kernels")
ENDIF()

include_directories(${TOP_SRC_DIR}/src)
#add_subdirectory(vigra_impex)
set(_COMMANDS enblend enfuse)

project(enblendenfuse)
#set(en_libraries vigra_impex ${common_libs})
set(en_libraries ${common_libs})

file(GLOB all_sources ${TOP_SRC_DIR}/src/*.cc)
foreach(_cmd ${_COMMANDS})
  list(REMOVE_ITEM all_sources "${TOP_SRC_DIR}/src/${_cmd}.cc")
endforeach(_cmd)

if (WIN32)
  file(GLOB winhelper_sources ${TOP_SRC_DIR}/src/win32helpers/*.c ${TOP_SRC_DIR}/src/win32helpers/*.h)
  list(APPEND all_sources ${winhelper_sources})
endif(WIN32)

file(GLOB layer_selection ${TOP_SRC_DIR}/src/layer_selection/*.cc ${TOP_SRC_DIR}/src/layer_selection/*.h)
list(APPEND all_sources ${layer_selection})
include_directories(${TOP_SRC_DIR}/src/layer_selection)
file(GLOB headers ${TOP_SRC_DIR}/src/*.h)
list(APPEND all_sources ${headers})
list(APPEND all_sources "${CMAKE_CURRENT_BINARY_DIR}/signature.h")
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(depends2)
message(STATUS "Boost_FOUND = ${Boost_FOUND}")
message(STATUS "OpenMP_CXX_FLAGS = ${OpenMP_CXX_FLAGS}")
foreach(_cmd ${_COMMANDS})
  set(${_cmd}_sources ${all_sources} "${TOP_SRC_DIR}/src/${_cmd}.cc")
  add_executable(${_cmd} ${${_cmd}_sources} )
  #message(STATUS "depends of ${_cmd} ${all_sources}")
  if(OpenMP_CXX_FLAGS AND NOT MSVC)
    set_target_properties(${_cmd} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
    message(STATUS "Adding PROPERTIES LINK_FLAGS to ${_cmd}")
  endif(OpenMP_CXX_FLAGS AND NOT MSVC)
  target_link_libraries(${_cmd} ${en_libraries})
  install(TARGETS ${_cmd} DESTINATION bin CONFIGURATIONS Release RelWithDebInfo MinSizeRel)
  #add_dependencies(${_cmd} vigra_impex)
  IF(NOT ${PREFER_SEPARATE_OPENCL_SOURCE})
    add_dependencies(${_cmd} cl_sources)
  ENDIF()
  list(APPEND depends2 ${_cmd})
endforeach(_cmd)
add_custom_target(enblendenfuse ALL DEPENDS ${depends2})

project(man)
# create enblend.1 and enfuse.1

find_program(HELP2MAN_EXE "help2man")
set(depends3)
foreach (_cmd ${_COMMANDS})
  if (NOT DOC OR ${HELP2MAN_EXE} MATCHES "-NOTFOUND")
    add_custom_command(
      OUTPUT "${_cmd}.1"
      COMMAND ${CMAKE_COMMAND} -E copy "${TOP_SRC_DIR}/src/${_cmd}.1" "${_cmd}.1"
      DEPENDS "${TOP_SRC_DIR}/src/${_cmd}.1"
      )
  else()
    add_custom_command(
      OUTPUT "${_cmd}.1"
      COMMAND "${HELP2MAN_EXE}" "--output=${_cmd}.1" "${CMAKE_BINARY_DIR}/bin/${_cmd}"
      DEPENDS "${CMAKE_BINARY_DIR}/bin/${_cmd}"
      )
  endif()
  list(APPEND depends3 "${CMAKE_CURRENT_BINARY_DIR}/${_cmd}.1")
endforeach(_cmd)
add_custom_target(man ALL DEPENDS ${depends3})
add_dependencies(man enblendenfuse)
