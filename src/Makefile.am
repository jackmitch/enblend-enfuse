SUBDIRS = vigra_impex
bin_PROGRAMS = enblend
enblend_SOURCES = assemble.cc blend.cc bounds.cc enblend.cc io.cc mask.cc nearest.cc pyramid.cc enblend.h
enblend_LDADD = vigra_impex/libvigra_impex.a -lpng -ljpeg -ltiff
enblend_CXXFLAGS = -g -O3 -Wall -I../include
