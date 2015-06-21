AUTOMAKE_OPTIONS = foreign no-installinfo

EXTRA_DIST = \
        README.txt VERSION \
        CMakeLists.txt CMakeModules ConfigureChecks.cmake config.h.cmake

#--if BUILD_DOC
#--  MAYBE_DOC = doc
#--endif

SUBDIRS = share src $(MAYBE_DOC)
DIST_SUBDIRS = share src

ACLOCAL_AMFLAGS = -I m4

# created by configure(1)
DISTCLEANFILES = config-h.texi
