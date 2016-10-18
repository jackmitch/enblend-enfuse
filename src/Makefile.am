AM_CPPFLAGS = $(EXTRACPPFLAGS)
AM_CFLAGS = $(EXTRACFLAGS)
AM_CXXFLAGS = $(EXTRACXXFLAGS)
AM_LDFLAGS = $(EXTRALDFLAGS)

SUBDIRS = dynamic_loader layer_selection win32helpers

bin_PROGRAMS = enblend enfuse

enblend_SOURCES = fillpolygon.hxx functoraccessor.hxx rect2d.hxx stride.hxx \
                  \
                  allocate.h \
                  anneal.h assemble.h blend.h bounds.h \
                  common.h enblend.h enblend.cc fixmath.h \
                  global.h graphcut.h \
                  maskcommon.h masktypedefs.h mask.h postoptimizer.h \
                  nearest.h numerictraits.h \
                  opencl.h opencl.cc opencl_anneal.h opencl_vigra.h \
                  openmp_def.h openmp_lock.h openmp_vigra.h \
                  path.h pyramid.h \
                  alternativepercentage.h alternativepercentage.cc \
                  error_message.h error_message.cc \
                  filenameparse.h filenameparse.cc \
                  filespec.h filespec.cc \
                  introspection.h introspection.cc \
                  mersenne.h mersenne.cc \
                  metadata.h metadata.cc \
                  parameter.h parameter.cc \
                  self_test.h self_test.cc \
                  tiff_message.h tiff_message.cc \
                  timer.h timer.cc \
                  minimizer.h minimizer.cc \
                  muopt.h optional_transitional.hpp
enblend_LDFLAGS = $(AM_LDFLAGS)
enblend_LDADD = layer_selection/liblayersel.a \
                $(GSL_LIBS) $(STATIC_LIBS) \
                @EXTRA_LIBS@
enblend_CXXFLAGS = -DENBLEND_SOURCE $(AM_CXXFLAGS) $(GSL_CFLAGS) \
                   -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 \
                   -I$(top_srcdir)/src/dynamic_loader \
                   -I$(top_srcdir)/src/layer_selection

enfuse_SOURCES = functoraccessor.hxx rect2d.hxx stride.hxx \
                 \
                 allocate.h \
                 assemble.h blend.h bounds.h common.h \
                 exposure_weight_base.h \
                 exposure_weight.h exposure_weight.cc \
                 enfuse.h enfuse.cc fixmath.h \
                 global.h mga.h numerictraits.h \
                 opencl.h opencl.cc opencl_vigra.h \
                 opencl_exposure_weight.h opencl_exposure_weight.cc \
                 openmp_def.h openmp_lock.h openmp_vigra.h \
                 pyramid.h \
                 alternativepercentage.h alternativepercentage.cc \
                 error_message.h error_message.cc \
                 filenameparse.h filenameparse.cc \
                 filespec.h filespec.cc \
                 introspection.h introspection.cc \
                 mersenne.h mersenne.cc \
                 metadata.h metadata.cc \
                 parameter.h parameter.cc \
                 self_test.h self_test.cc \
                 tiff_message.h tiff_message.cc \
                 timer.h timer.cc \
                 minimizer.h minimizer.cc \
                 muopt.h optional_transitional.hpp
enfuse_LDFLAGS = $(AM_LDFLAGS)
enfuse_LDADD = dynamic_loader/libdynamic_loader.a \
               layer_selection/liblayersel.a \
               $(GSL_LIBS) $(STATIC_LIBS) \
               @EXTRA_LIBS@
enfuse_CXXFLAGS = -DENFUSE_SOURCE $(AM_CXXFLAGS) $(GSL_CFLAGS) \
                  -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 \
                  -I$(top_srcdir)/src/dynamic_loader \
                  -I$(top_srcdir)/src/layer_selection

EXTRA_DIST = embrace calculate_state_probabilities.cl distance_transform_fh.cl \
             enblend.1 enfuse.1 \
             gen_sig DefaultSig.pm Sig.pm \
             CMakeLists.txt

DISTCLEANFILES = enblend.1 enfuse.1


# Generated sources

BUILT_SOURCES = signature.h calculate_state_probabilities.icl distance_transform_fh.icl
CLEANFILES = signature.h calculate_state_probabilities.icl distance_transform_fh.icl

%.icl: %.cl
	$(AM_V_GEN)$(srcdir)/embrace --format=c++ --label=$(basename $(notdir $<))_source_code $< > $@

calculate_state_probabilities.icl: calculate_state_probabilities.cl
distance_transform_fh.icl: distance_transform_fh.cl

signature.h: $(srcdir)/gen_sig $(srcdir)/DefaultSig.pm $(srcdir)/Sig.pm
	@ $(PERL) -I$(srcdir) $< --extra=$(VERSION) > $@


# Documentation

man_MANS = enblend.1 enfuse.1
#dist_man_MANS = enblend.1 enfuse.1

.PHONY: man
man: enblend.1 enfuse.1


enblend.1: enblend
	$(AM_V_GEN)$(HELP2MAN) --no-info --output=$@ ./enblend

enfuse.1: enfuse
	$(AM_V_GEN)$(HELP2MAN) --no-info --output=$@ ./enfuse
