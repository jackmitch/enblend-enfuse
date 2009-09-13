# This file is part of enblend/enfuse.
# Licence details can be found in the file COPYING.
#
# Copyright (c) 2009, Kornel Benko <Kornel.Benko@berlin.de>
#                   , Ryan Sleevi <ryan+hugin@sleevi.com>
#                   , Harry van der Wolf <hvdwolf@gmail.com>
#

include_directories(${TOP_SRC_DIR}/src)
add_subdirectory(vigra_impex)

project(enblendenfuse)
file(GLOB all_sources ${TOP_SRC_DIR}/src/*.cc)
list(REMOVE_ITEM all_sources ${TOP_SRC_DIR}/src/enblend.cc ${TOP_SRC_DIR}/src/enfuse.cc)

if (WIN32)
  file(GLOB winhelper_sources ${TOP_SRC_DIR}/src/win32helpers/*.c)
  list(APPEND all_sources ${winhelper_sources})
endif(WIN32)


set(enblend_sources ${all_sources} ${TOP_SRC_DIR}/src/enblend.cc)
set(enfuse_sources ${all_sources} ${TOP_SRC_DIR}/src/enfuse.cc)
add_executable(enblend ${enblend_sources})
add_executable(enfuse ${enfuse_sources})

if (${OpenMP_CXX_FLAGS})
  set_target_properties(enblend enfuse PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
endif(${OpenMP_CXX_FLAGS})

set(en_libraries vigra_impex ${common_libs})
target_link_libraries(enblend ${en_libraries})
target_link_libraries(enfuse ${en_libraries})

install(TARGETS enblend enfuse DESTINATION bin CONFIGURATIONS Release RelWithDebInfo MinSizeRel)
add_custom_target(enblendenfuse ALL DEPENDS enblend enfuse)

add_dependencies(enblend vigra_impex)
add_dependencies(enfuse vigra_impex)

project(man)
# create enblend.1 and enfuse.1

find_program(HELP2MAN_EXE "help2man")
foreach (_cmd "enblend" "enfuse")
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
  add_custom_target(man ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${_cmd}.1)
endforeach(_cmd)

add_dependencies(man enblendenfuse)
