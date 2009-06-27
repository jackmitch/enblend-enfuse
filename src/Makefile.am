SUBDIRS = vigra_impex win32helpers

bin_PROGRAMS = enblend enfuse

enblend_SOURCES = anneal.h assemble.h blend.h bounds.h \
                  common.h enblend.h enblend.cc fixmath.h \
                  float_cast.h global.h gpu.cc gpu.h mask.h \
                  nearest.h numerictraits.h path.h pyramid.h \
                  filenameparse.h filenameparse.cc
enblend_LDFLAGS =
enblend_LDADD = vigra_impex/libvigra_impex.a @EXTRA_LIBS@
enblend_CXXFLAGS = -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -I${top_srcdir}/include

enfuse_SOURCES = assemble.h blend.h bounds.h common.h \
                 enfuse.h enfuse.cc fixmath.h float_cast.h \
                 global.h mga.h numerictraits.h pyramid.h \
                 filenameparse.h filenameparse.cc
enfuse_LDFLAGS =
enfuse_LDADD = vigra_impex/libvigra_impex.a @EXTRA_LIBS@
enfuse_CXXFLAGS = -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -I${top_srcdir}/include

EXTRA_DIST = enblend.vcproj enblend_msvc2003.vcproj enblend_msvc2008.vcproj \
             enfuse.vcproj enfuse_msvc2003.vcproj enfuse_msvc2008.vcproj \
             enblend.1 enfuse.1

DISTCLEANFILES = enblend.1 enfuse.1

# Documentation

man_MANS = enblend.1 enfuse.1
#dist_man_MANS = enblend.1 enfuse.1

.PHONY: man
man: enblend.1 enfuse.1


enblend.1: enblend
	$(HELP2MAN) --output=$@ ./enblend

enfuse.1: enfuse
	$(HELP2MAN) --output=$@ ./enfuse
