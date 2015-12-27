EXTRA_DIST = \
        README.txt VERSION \
        CMakeLists.txt CMakeModules ConfigureChecks.cmake config.h.cmake

if BUILD_DOC
  MAYBE_DOC = doc
endif

SUBDIRS = src $(MAYBE_DOC)
DIST_SUBDIRS = src doc

ACLOCAL_AMFLAGS = -I m4
