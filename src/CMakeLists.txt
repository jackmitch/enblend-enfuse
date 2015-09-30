# This file is part of enblend/enfuse.
# Licence details can be found in the file COPYING.
#
# Copyright (c) 2009-2015, Kornel Benko <Kornel.Benko@berlin.de>
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
add_custom_target(signature DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/signature.h")

IF(ENABLE_OPENCL)
  # prepare all OpenCL sources
  SET(OPENCL_SOURCES)
  FILE(GLOB opencl_sources ${TOP_SRC_DIR}/src/*.cl)
  IF(PREFER_SEPARATE_OPENCL_SOURCE)
    INSTALL(FILES ${opencl_sources} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/enblend/kernels)
  ELSE()
    FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/kernels")
    FOREACH(_cl_source ${opencl_sources})
      GET_FILENAME_COMPONENT(_cl_filename "${_cl_source}" NAME)

      STRING(REGEX REPLACE "[.]cl$" ".icl" _icl_source "${_cl_filename}")
      SET(_icl_source "${CMAKE_CURRENT_BINARY_DIR}/kernels/${_icl_source}")
      STRING(REGEX REPLACE "[.]cl$" "_source_code" _cl_label "${_cl_filename}")

      ADD_CUSTOM_COMMAND(
        OUTPUT "${_icl_source}"
        COMMAND ${PERL_EXECUTABLE} "${TOP_SRC_DIR}/src/embrace"
          "--format=c++" "--label=${_cl_label}" "${_cl_source}" ">" "${_icl_source}"
        DEPENDS "${_cl_source}"
        COMMENT "Generating ${_icl_source}"
      )
      SET_SOURCE_FILES_PROPERTIES("${_icl_source}" GENERATED)
      LIST(APPEND OPENCL_SOURCES "${_icl_source}")
    ENDFOREACH()

    ADD_CUSTOM_TARGET(cl_sources ALL DEPENDS "${OPENCL_SOURCES}")
    INCLUDE_DIRECTORIES("${CMAKE_CURRENT_BINARY_DIR}/kernels")
  ENDIF()
ENDIF()

include_directories(${TOP_SRC_DIR}/src)
set(ENBLEND_SOURCES 
    fillpolygon.hxx functoraccessor.hxx rect2d.hxx stride.hxx 
    anneal.h assemble.h blend.h bounds.h
    common.h enblend.h enblend.cc fixmath.h
    global.h graphcut.h
    maskcommon.h masktypedefs.h mask.h postoptimizer.h
    nearest.h numerictraits.h
    opencl.h opencl.cc opencl_vigra.h
    openmp_def.h openmp_lock.h openmp_vigra.h
    path.h pyramid.h
    alternativepercentage.h alternativepercentage.cc
    error_message.h error_message.cc
    filenameparse.h filenameparse.cc
    filespec.h filespec.cc
    introspection.h introspection.cc
    mersenne.h mersenne.cc
    parameter.h parameter.cc
    self_test.h self_test.cc
    tiff_message.h tiff_message.cc
    timer.h timer.cc
    minimizer.h minimizer.cc
    muopt.h
)
set(ENFUSE_SOURCES 
    functoraccessor.hxx rect2d.hxx stride.hxx
    assemble.h blend.h bounds.h common.h
    exposure_weight_base.h
    exposure_weight.h exposure_weight.cc
    enfuse.h enfuse.cc fixmath.h
    global.h mga.h numerictraits.h
    opencl.h opencl.cc opencl_vigra.h
    openmp_def.h openmp_lock.h openmp_vigra.h
    pyramid.h
    alternativepercentage.h alternativepercentage.cc
    error_message.h error_message.cc
    filenameparse.h filenameparse.cc
    filespec.h filespec.cc
    introspection.h introspection.cc
    mersenne.h mersenne.cc
    parameter.h parameter.cc
    self_test.h self_test.cc
    tiff_message.h tiff_message.cc
    timer.h timer.cc
    minimizer.h minimizer.cc
    muopt.h
)

if (WIN32)
  add_library(getopt STATIC win32helpers/getopt.c win32helpers/getopt.h)
  set(additional_libs getopt)
else()
  set(additional_libs)
endif()

set(LAYER_SELECTION_FILES
    layer_selection/info.h layer_selection/info.cc
    layer_selection/layer_selection.h layer_selection/layer_selection.cc
    layer_selection/selector.h layer_selection/selector.cc
)
add_library(layer_selection STATIC ${LAYER_SELECTION_FILES})
include_directories(${TOP_SRC_DIR}/src/layer_selection)

set(DYNAMIC_LOADER_FILES
    dynamic_loader/dynamic_loader.h
    dynamic_loader/dynamic_loader.cc
)
add_library(dynamic_loader STATIC ${DYNAMIC_LOADER_FILES})
include_directories(${TOP_SRC_DIR}/src/dynamic_loader)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

message(STATUS "Boost_FOUND = ${Boost_FOUND}")
message(STATUS "OpenMP_CXX_FLAGS = ${OpenMP_CXX_FLAGS}")

add_executable(enblend ${ENBLEND_SOURCES})
add_executable(enfuse ${ENFUSE_SOURCES})

add_dependencies(enblend signature)
add_dependencies(enfuse signature)
target_compile_definitions(enblend PUBLIC "-DENBLEND_SOURCE")
target_compile_definitions(enfuse PUBLIC "-DENFUSE_SOURCE")
if(OpenMP_CXX_FLAGS AND NOT MSVC)
    set_target_properties(enblend PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
    set_target_properties(enfuse PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
    message(STATUS "Adding PROPERTIES LINK_FLAGS to enblend and enfuse")
endif(OpenMP_CXX_FLAGS AND NOT MSVC)
target_link_libraries(enblend ${common_libs} ${additional_libs} layer_selection)
target_link_libraries(enfuse ${common_libs} ${additional_libs} layer_selection dynamic_loader)
install(TARGETS enblend enfuse DESTINATION bin CONFIGURATIONS Release RelWithDebInfo MinSizeRel)

IF(ENABLE_OPENCL AND NOT ${PREFER_SEPARATE_OPENCL_SOURCE})
    add_dependencies(enblend cl_sources)
    add_dependencies(enfuse cl_sources)
ENDIF()

if(NOT WIN32)
    # create enblend.1 and enfuse.1
    if(NOT MANDIR AND NOT $ENV{MANDIR} STREQUAL "")
        set(MANDIR "$ENV{MANDIR}")
    elseif(CMAKE_SYSTEM_NAME STREQUAL "FreeBSD" OR CMAKE_SYSTEM_NAME STREQUAL "NetBSD")
        set(MANDIR "man")
    else(NOT MANDIR AND NOT $ENV{MANDIR} STREQUAL "")
        set(MANDIR "share/man")
    endif()
    
    find_program(HELP2MAN_EXE "help2man")
    if (NOT DOC OR ${HELP2MAN_EXE} MATCHES "-NOTFOUND")
        INSTALL(FILES enblend.1 enfuse.1 DESTINATION ${MANDIR}/man1)
    else()
        project(man)
        add_custom_command(
          OUTPUT "${CMAKE_BINARY_DIR}/src/enblend.1"
          COMMAND "${HELP2MAN_EXE}" "--output=enblend.1" "${CMAKE_BINARY_DIR}/bin/enblend"
          DEPENDS "${CMAKE_BINARY_DIR}/bin/enblend"
          )
        add_custom_command(
          OUTPUT "${CMAKE_BINARY_DIR}/src/enfuse.1"
          COMMAND "${HELP2MAN_EXE}" "--output=enfuse.1" "${CMAKE_BINARY_DIR}/bin/enfuse"
          DEPENDS "${CMAKE_BINARY_DIR}/bin/enfuse"
          )
        add_custom_target(man ALL DEPENDS enblend.1 enfuse.1)
        add_dependencies(man enblend enfuse)
        INSTALL(FILES ${CMAKE_BINARY_DIR}/src/enblend.1 ${CMAKE_BINARY_DIR}/src/enfuse.1 DESTINATION ${MANDIR}/man1)
    endif()
endif()
