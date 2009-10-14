EXTRA_DIST = VIGRA_LICENSE README.txt \
             enblend_msvc2008.sln

if BUILD_DOC
  MAYBE_DOC = doc
endif

SUBDIRS = include src $(MAYBE_DOC)
DIST_SUBDIRS = include src doc

ACLOCAL_AMFLAGS = -I m4

# created by configure(1)
DISTCLEANFILES = config-h.texi
