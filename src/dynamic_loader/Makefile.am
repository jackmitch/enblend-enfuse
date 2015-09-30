AM_CPPFLAGS = $(EXTRACPPFLAGS)
AM_CFLAGS = $(EXTRACFLAGS)
AM_CXXFLAGS = $(EXTRACXXFLAGS)
AM_LDFLAGS = $(EXTRALDFLAGS)

noinst_LIBRARIES = libdynamic_loader.a

libdynamic_loader_a_SOURCES = dynamic_loader.h dynamic_loader.cc \
                              dynamic_loader_implementation.h dynamic_loader_implementation.cc \
                              gmodule_implementation.h gmodule_implementation.cc \
                              null_implementation.h \
                              sunny_implementation.h sunny_implementation.cc \
                              win32_implementation.h win32_implementation.cc

libdynamic_loader_a_CXXFLAGS = $(AM_CXXFLAGS) \
                         -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 \
                         -I$(top_srcdir)/include -I$(top_srcdir)/src -I$(top_srcdir)/src/dynamic_loader
