EXTRA_DIST = enblend.vcproj enblend_msvc2003.vcproj enblend_msvc2008.vcproj enfuse_msvc2003.vcproj enfuse_msvc2003.vcproj
SUBDIRS = vigra_impex win32helpers
bin_PROGRAMS = enblend enfuse
enblend_SOURCES = anneal.h assemble.h blend.h bounds.h common.h enblend.h enblend.cc fixmath.h float_cast.h gpu.cc gpu.h mask.h nearest.h numerictraits.h path.h pyramid.h
enblend_LDFLAGS = -DENBLEND_CACHE_IMAGES -ffast-math
#enblend_LDADD = vigra_impex/libvigra_impex.a -lpng -ljpeg -ltiff
enblend_LDADD = vigra_impex/libvigra_impex.a
enblend_CXXFLAGS = -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -DENBLEND_CACHE_IMAGES -I${top_srcdir}/include
enfuse_SOURCES = assemble.h blend.h bounds.h common.h enfuse.h enfuse.cc fixmath.h float_cast.h numerictraits.h pyramid.h
enfuse_LDFLAGS = -DENBLEND_CACHE_IMAGES -ffast-math
enfuse_LDADD = vigra_impex/libvigra_impex.a
enfuse_CXXFLAGS = -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -DENBLEND_CACHE_IMAGES -I${top_srcdir}/include
