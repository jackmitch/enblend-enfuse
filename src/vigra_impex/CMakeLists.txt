# This file is part of enblend/enfuse.
# Licence details can be found in the file COPYING.
#
# Copyright (c) 2009-2010, Kornel Benko <Kornel.Benko@berlin.de>
#

project(vigra_impex)
file(GLOB vigra_impex_cxx_sources ${TOP_SRC_DIR}/src/vigra_impex/*.cxx ${TOP_SRC_DIR}/src/vigra_impex/*.c)
file(GLOB vigra_impex_hxx_headers ${TOP_SRC_DIR}/src/vigra_impex/*.hxx ${TOP_SRC_DIR}/src/vigra_impex/*.h)

#add_definitions(${EXR_CFLAGS})
include_directories(${TOP_SRC_DIR}/src/vigra_impex)
IF(OPENEXR_FOUND)
include_directories(${OPENEXR_INCLUDE_DIR})
ENDIF()
add_library(vigra_impex ${library_type} ${vigra_impex_cxx_sources})
target_link_libraries(vigra_impex ${JPEG_LIBRARIES} ${OPENEXR_LIBRARIES} ${PNG_LIBRARIES} ${TIFF_LIBRARIES})
