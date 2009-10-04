AM_CPPFLAGS = $(EXTRACPPFLAGS)
AM_CFLAGS = $(EXTRACFLAGS)
AM_CXXFLAGS = $(EXTRACXXFLAGS)
AM_LDFLAGS = $(EXTRALDFLAGS)

SUBDIRS = vigra_impex win32helpers

bin_PROGRAMS = enblend enfuse

enblend_SOURCES = anneal.h assemble.h blend.h bounds.h \
                  common.h enblend.h enblend.cc fixmath.h \
                  float_cast.h global.h gpu.cc gpu.h mask.h \
                  nearest.h numerictraits.h openmp.h path.h pyramid.h \
                  filenameparse.h filenameparse.cc
enblend_LDFLAGS = $(AM_LDFLAGS) $(OPENGL_CFLAGS)
enblend_LDADD = vigra_impex/libvigra_impex.a $(OPENGL_LIBS) @EXTRA_LIBS@
enblend_CXXFLAGS = $(AM_CXXFLAGS) $(OPENGL_CFLAGS) \
                   -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 \
                   -I${top_srcdir}/include

enfuse_SOURCES = assemble.h blend.h bounds.h common.h \
                 enfuse.h enfuse.cc fixmath.h float_cast.h \
                 global.h mga.h numerictraits.h openmp.h pyramid.h \
                 filenameparse.h filenameparse.cc
enfuse_LDFLAGS = $(AM_LDFLAGS)
enfuse_LDADD = vigra_impex/libvigra_impex.a @EXTRA_LIBS@
enfuse_CXXFLAGS = $(AM_CXXFLAGS) \
                  -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 \
                  -I${top_srcdir}/include

EXTRA_DIST = enblend.vcproj enblend_msvc2003.vcproj enblend_msvc2008.vcproj \
             enfuse.vcproj enfuse_msvc2003.vcproj enfuse_msvc2008.vcproj \
             enblend.1 enfuse.1 \
             gen_sig DefaultSig.pm Sig.pm

DISTCLEANFILES = enblend.1 enfuse.1


# Generated sources

BUILT_SOURCES = signature.h
CLEANFILES = signature.h

signature.h: $(srcdir)/gen_sig $(srcdir)/DefaultSig.pm $(srcdir)/Sig.pm
	@ $(PERL) -I$(srcdir) $< --extra=$(VERSION) > $@


# Documentation

man_MANS = enblend.1 enfuse.1
#dist_man_MANS = enblend.1 enfuse.1

.PHONY: man
man: enblend.1 enfuse.1


enblend.1: enblend
	$(HELP2MAN) --output=$@ ./enblend

enfuse.1: enfuse
	$(HELP2MAN) --output=$@ ./enfuse
