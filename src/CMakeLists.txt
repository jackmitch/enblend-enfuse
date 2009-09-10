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
message(STATUS "all_sources = " ${all_sources})
list(REMOVE_ITEM all_sources ${TOP_SRC_DIR}/src/enblend.cc ${TOP_SRC_DIR}/src/enfuse.cc)
if (WIN32)
    list(APPEND all_sources ${TOP_SRC_DIR}/src/win32helpers/getopt.c ${TOP_SRC_DIR}/src/win32helpers/getopt_long.c)
endif()

message(STATUS "all_sources(after remove) = " ${all_sources})
set(enblend_sources ${all_sources} ${TOP_SRC_DIR}/src/enblend.cc)
set(enfuse_sources ${all_sources} ${TOP_SRC_DIR}/src/enfuse.cc)
add_executable(enblend ${enblend_sources})
add_executable(enfuse ${enfuse_sources})
message(STATUS "OPENMP_FOUND = ${OPENMP_FOUND}")
if (OPENMP_FOUND)
  set_target_properties(enblend enfuse PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
  message(STATUS "OpenMP_CXX_FLAGS = ${OpenMP_CXX_FLAGS}")
endif()

set(en_libraries vigra_impex ${common_libs})
target_link_libraries(enblend ${en_libraries})
target_link_libraries(enfuse ${en_libraries})

