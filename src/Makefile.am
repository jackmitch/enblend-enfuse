SUBDIRS = vigra_impex
bin_PROGRAMS = enblend
enblend_SOURCES = alpha.cc assemble.h common.h enblend.h enblend.cc
enblend_LDADD = vigra_impex/libvigra_impex.a -lpng -ljpeg -ltiff
enblend_CXXFLAGS = -g -O3 -Wall -I../include
