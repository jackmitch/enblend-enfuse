# This file is part of enblend/enfuse.
# Licence details can be found in the file COPYING.
#
# Copyright (c) 2009, Kornel Benko <Kornel.Benko@berlin.de>
#                   , Ryan Sleevi <ryan+hugin@sleevi.com>
#                   , Harry van der Wolf <hvdwolf@gmail.com>
#

include_directories(${TOP_SRC_DIR}/src)
add_subdirectory(vigra_impex)
set(_COMMANDS enblend enfuse)

project(enblendenfuse)
set(en_libraries vigra_impex ${common_libs})

file(GLOB all_sources ${TOP_SRC_DIR}/src/*.cc)
foreach(_cmd ${_COMMANDS})
  list(REMOVE_ITEM all_sources "${TOP_SRC_DIR}/src/${_cmd}.cc")
endforeach(_cmd)

if (WIN32)
  file(GLOB winhelper_sources ${TOP_SRC_DIR}/src/win32helpers/*.c)
  list(APPEND all_sources ${winhelper_sources})
endif(WIN32)

foreach(_cmd ${_COMMANDS})
  set(${_cmd}_sources ${all_sources} "${TOP_SRC_DIR}/src/${_cmd}.cc")
  add_executable(${_cmd} ${${_cmd}_sources})
  if (${OpenMP_CXX_FLAGS})
    set_target_properties(${_cmd} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif(${OpenMP_CXX_FLAGS})
  target_link_libraries(${_cmd} ${en_libraries})
  install(TARGETS ${_cmd} DESTINATION bin CONFIGURATIONS Release RelWithDebInfo MinSizeRel)
  add_dependencies(${_cmd} vigra_impex)

  add_custom_target(enblendenfuse ALL DEPENDS ${_cmd})
endforeach(_cmd)

project(man)
# create enblend.1 and enfuse.1

find_program(HELP2MAN_EXE "help2man")
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
  add_custom_target(man ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${_cmd}.1)
endforeach(_cmd)

add_dependencies(man enblendenfuse)
