# This file is part of enblend/enfuse.
# Licence details can be found in the file COPYING.
#
# Copyright (c) 2009, Kornel Benko <Kornel.Benko@berlin.de>
#                   , Ryan Sleevi <ryan+hugin@sleevi.com>
#                   , Harry van der Wolf <hvdwolf@gmail.com>
#

include_directories(${TOP_SRC_DIR}/src)
add_subdirectory(vigra_impex)

file(GLOB all_sources ${TOP_SRC_DIR}/src/*.cc)
list(REMOVE_ITEM all_sources ${TOP_SRC_DIR}/src/enblend.cc ${TOP_SRC_DIR}/src/enfuse.cc)

if (WIN32)
  list(APPEND all_sources ${TOP_SRC_DIR}/src/win32helpers/getopt.c ${TOP_SRC_DIR}/src/win32helpers/getopt_long.c)
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

# create enblend.1 and enfuse.1

find_program(HELP2MAN_EXE "help2man")
foreach (_cmd "enblend" "enfuse")
  set(_binary "${CMAKE_BINARY_DIR}/bin/${_cmd}")
  if (${HELP2MAN_EXE} MATCHES "-NOTFOUND")
    add_custom_command(
      OUTPUT "${_cmd}.1"
      COMMAND ${CMAKE_COMMAND} -E copy "${TOP_SRC_DIR}/src/${_cmd}.1" "${_cmd}.1"
      DEPENDS "${TOP_SRC_DIR}/src/${_cmd}.1"
      )
  else()
    add_custom_command(
      OUTPUT "${_cmd}.1"
      COMMAND "${HELP2MAN_EXE}" "--output=${_cmd}.1" "${_binary}"
      DEPENDS "${_binary}"
      )
  endif()
  install(FILES "${_cmd}.1" DESTINATION "share/man/man1")
  ADD_CUSTOM_TARGET("enblend" ALL DEPENDS "${_cmd}.1")
endforeach(_cmd)

