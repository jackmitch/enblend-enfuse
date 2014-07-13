AM_CPPFLAGS = $(EXTRACPPFLAGS)
AM_CFLAGS = $(EXTRACFLAGS)
AM_CXXFLAGS = $(EXTRACXXFLAGS)
AM_LDFLAGS = $(EXTRALDFLAGS)

SUBDIRS = layer_selection win32helpers

bin_PROGRAMS = enblend enfuse

enblend_SOURCES = fillpolygon.hxx functoraccessor.hxx rect2d.hxx stride.hxx \
                  \
                  anneal.h assemble.h blend.h bounds.h \
                  common.h enblend.h enblend.cc fixmath.h \
                  global.h graphcut.h \
                  maskcommon.h masktypedefs.h mask.h postoptimizer.h \
                  nearest.h numerictraits.h \
                  opencl.h opencl.cc opencl_vigra.h \
                  openmp_def.h openmp_lock.h openmp_vigra.h \
                  path.h pyramid.h \
                  error_message.h error_message.cc \
                  filenameparse.h filenameparse.cc \
                  filespec.h filespec.cc \
                  self_test.h self_test.cc \
                  tiff_message.h tiff_message.cc \
                  timer.h timer.cc \
                  minimizer.h muopt.h
enblend_LDFLAGS = $(AM_LDFLAGS)
enblend_LDADD = layer_selection/liblayersel.a \
                $(GSL_LIBS) @EXTRA_LIBS@
enblend_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CFLAGS) \
                   -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 \
                   -I${top_srcdir}/src/layer_selection

enfuse_SOURCES = functoraccessor.hxx rect2d.hxx stride.hxx \
                 \
                 assemble.h blend.h bounds.h common.h \
                 exposure_weight_base.h exposure_weight.h dynamic_loader.h \
                 enfuse.h enfuse.cc fixmath.h \
                 global.h mga.h numerictraits.h \
                 opencl.h opencl.cc opencl_vigra.h \
                 openmp_def.h openmp_lock.h openmp_vigra.h \
                 pyramid.h \
                 error_message.h error_message.cc \
                 filenameparse.h filenameparse.cc \
                 filespec.h filespec.cc \
                 self_test.h self_test.cc \
                 tiff_message.h tiff_message.cc \
                 timer.h timer.cc \
                 minimizer.h muopt.h
enfuse_LDFLAGS = $(AM_LDFLAGS)
enfuse_LDADD = layer_selection/liblayersel.a \
               $(GSL_LIBS) @EXTRA_LIBS@
enfuse_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CFLAGS) \
                  -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 \
                  -I${top_srcdir}/src/layer_selection

EXTRA_DIST = embrace distance_transform_fh.cl \
             enblend.1 enfuse.1 \
             gen_sig DefaultSig.pm Sig.pm \
             CMakeLists.txt

DISTCLEANFILES = enblend.1 enfuse.1


# Generated sources

BUILT_SOURCES = signature.h distance_transform_fh.icl
CLEANFILES = signature.h distance_transform_fh.icl

signature.h: $(srcdir)/gen_sig $(srcdir)/DefaultSig.pm $(srcdir)/Sig.pm
	@ $(PERL) -I$(srcdir) $< --extra=$(VERSION) > $@

distance_transform_fh.icl: distance_transform_fh.cl
	$(srcdir)/embrace --format=c++ --label=distance_transform_fh_source_code $< > $@


# Documentation

man_MANS = enblend.1 enfuse.1
#dist_man_MANS = enblend.1 enfuse.1

.PHONY: man
man: enblend.1 enfuse.1


enblend.1: enblend
	$(HELP2MAN) --output=$@ ./enblend

enfuse.1: enfuse
	$(HELP2MAN) --output=$@ ./enfuse
