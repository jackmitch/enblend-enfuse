AM_CPPFLAGS = $(EXTRACPPFLAGS)
AM_CFLAGS = $(EXTRACFLAGS)
AM_CXXFLAGS = $(EXTRACXXFLAGS)
AM_LDFLAGS = $(EXTRALDFLAGS)

noinst_LIBRARIES = libvigra_impex.a

libvigra_impex_a_SOURCES = auto_file.hxx bmp.cxx bmp.hxx byteorder.cxx \
                           byteorder.hxx codecmanager.cxx codecmanager.hxx error.hxx \
                           gif.cxx gif.hxx hdr.cxx hdr.hxx \
                           iccjpeg.c iccjpeg.h imageinfo.cxx jpeg.cxx \
                           jpeg.hxx png.cxx png.hxx pnm.cxx \
                           pnm.hxx rgbe.c rgbe.h sun.cxx \
                           sun.hxx tiff.cxx tiff.hxx viff.cxx \
                           viff.hxx void_vector.cxx void_vector.hxx exr.cxx \
                           exr.hxx

libvigra_impex_a_CXXFLAGS = $(AM_CXXFLAGS) \
                            -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 \
                            -I${top_srcdir}/include
