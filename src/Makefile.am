SUBDIRS = vigra_impex
bin_PROGRAMS = enblend
enblend_SOURCES = assemble.h blend.h bounds.h common.h enblend.h enblend.cc fixmath.h mask.h nearest.h pyramid.h
#enblend_LDADD = vigra_impex/libvigra_impex.a -lpng -ljpeg -ltiff
enblend_LDFLAGS = -DENBLEND_CACHE_IMAGES -ffast-math
#enblend_LDADD = vigra_impex/libvigra_impex.a -ltiff
enblend_LDADD = vigra_impex/libvigra_impex.a /home/mihal/src/tiff-v3.6.1/libtiff/libtiff.a
enblend_CXXFLAGS = -g -O3 -ffast-math -Wall -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -DENBLEND_CACHE_IMAGES -I../include
