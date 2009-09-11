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

