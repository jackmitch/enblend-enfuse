EXTRA_DIST = enblend.vcproj
SUBDIRS = vigra_impex win32helpers
bin_PROGRAMS = enblend
enblend_SOURCES = anneal.h assemble.h blend.h bounds.h common.h enblend.h enblend.cc fixmath.h mask.h nearest.h pyramid.h
enblend_LDFLAGS = -DENBLEND_CACHE_IMAGES -ffast-math
#enblend_LDADD = vigra_impex/libvigra_impex.a -lpng -ljpeg -ltiff
enblend_LDADD = vigra_impex/libvigra_impex.a -ltiff
enblend_CXXFLAGS = -g -O3 -ffast-math -Wall -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -DENBLEND_CACHE_IMAGES -I../include
