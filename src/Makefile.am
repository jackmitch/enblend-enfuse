SUBDIRS = vigra_impex
bin_PROGRAMS = enblend
enblend_SOURCES = assemble.h bounds.h common.h enblend.h enblend.cc fixmath.h mask.h nearest.h pyramid.h
enblend_LDADD = vigra_impex/libvigra_impex.a -lpng -ljpeg -ltiff
enblend_CXXFLAGS = -g -O3 -Wall -I../include
