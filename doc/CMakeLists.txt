# This file is part of enblend.
# Licence details can be found in the file COPYING.
#
# Copyright (c) 2009, Kornel Benko <Kornel.Benko@berlin.de>
#

project(doc)

set(_docs)
set(_enblend_TEXINFOS)
foreach(_e auxmac.texi auxmac.tex fdl.texi mask-template-characters.texi helpful-programs.texi
    tuning-memory-usage.texi understanding-masks.texi workflow.texi authors.texi
    photographic-workflow.fig default.css)
  list(APPEND _enblend_TEXINFOS "${TOP_SRC_DIR}/doc/${_e}")
endforeach(_e)

set(_enfuse_TEXINFOS)
foreach(_e auxmac.texi auxmac.tex fdl.texi mask-template-characters.texi helpful-programs.texi
    tuning-memory-usage.texi understanding-masks.texi workflow.texi authors.texi
    config-edge.gp config.gp entropy-cutoff.gp.in entropy.gp.in gaussian.gp.in
    laplacian-of-gaussian.gp.in sharp-edge.gp.in sharp-edge.data
    smooth-edge.gp.in smooth-edge.data local-analysis-window.fig
    photographic-workflow.fig focus-stack-decision-tree.fig default.css)
  list(APPEND _enfuse_TEXINFOS "${TOP_SRC_DIR}/doc/${_e}")
endforeach(_e)

# First search for makeinfo command and
# generate info and pdf for enblend + enfuse
find_program(MAKEINFO_EXE "makeinfo")
find_program(TEXI2DVI_EXE "texi2dvi")
foreach (_cmd "enblend" "enfuse")
  set(_created "${CMAKE_CURRENT_BINARY_DIR}/${_cmd}")
  if (${MAKEINFO_EXE} MATCHES "-NOTFOUND" OR ${TEXI2DVI_EXE} MATCHES "-NOTFOUND")
    # there is no makeinfo or no texi2dvi command,
    # this implies insufficient latex installation.
    # Therefore simply copy the data from source tree.
    foreach(_ext "info" "pdf")
      SET_SOURCE_FILES_PROPERTIES("${_created}.${_ext}" GENERATED)
      add_custom_command(
	OUTPUT "${_created}.${_ext}"
	COMMAND ${CMAKE_COMMAND} -E copy "${TOP_SRC_DIR}/doc/${_cmd}.${_ext}" "${_created}.${_ext}"
	DEPENDS "${TOP_SRC_DIR}/doc/${_cmd}.${_ext}"
	)
      list(APPEND _docs "${_created}.${_ext}")
      install(FILES "${_created}.${_ext}" DESTINATION "doc")
    endforeach()
  else()
    # here makeinfo and texi2dvi exist
    add_custom_command(
      OUTPUT "${_created}.info"
      COMMAND ${MAKEINFO_EXE} "--css-include=${TOP_SRC_DIR}/doc/default.css"
      "-I" "${TOP_SRC_DIR}/doc"
      -o "${_created}.info" "${TOP_SRC_DIR}/doc/${_cmd}.texi"
      DEPENDS "${TOP_SRC_DIR}/doc/${_cmd}.texi" ${_${_cmd}_TEXINFOS}
      )
    message(STATUS "create ${_created}.pdf")
    add_custom_command(
      OUTPUT "${_created}.pdf"
      COMMAND ${TEXI2DVI_EXE} "--pdf" "--batch" "${TOP_SRC_DIR}/doc/${_cmd}.texi"
      "${TOP_SRC_DIR}/doc/${_cmd}.texi"
      DEPENDS "${TOP_SRC_DIR}/doc/${_cmd}.texi" ${_${_cmd}_TEXINFOS}
      )
    install(FILES "${_created}.info" "${_created}.pdf" DESTINATION "doc")
  endif()
  install(FILES "${_created}.info" "${_created}.pdf" DESTINATION "doc")
  list(APPEND _docs "${_created}.info" "${_created}.pdf")
endforeach(_cmd)

# search for gnuplot and
# generate images (.txt, .pdf .eps and .png)
find_program(GNUPLOT_EXE "gnuplot")
if (${GNUPLOT_EXE} MATCHES "-NOTFOUND")
  # there is no gnuplot command, so simply copy the needed files
  file(GLOB _copy_files RELATIVE "${TOP_SRC_DIR}/doc"
    "${TOP_SRC_DIR}/doc/*.txt"
    "${TOP_SRC_DIR}/doc/*.eps"
    "${TOP_SRC_DIR}/doc/*.png"
    "${TOP_SRC_DIR}/doc/*.pdf")
  foreach(_cf ${_copy_files})
    add_custom_command(
      OUTPUT ${_cf}
      COMMAND ${CMAKE_COMMAND} -E copy "${TOP_SRC_DIR}/doc/${_cf}" "${_cf}"
      DEPENDS "${TOP_SRC_DIR}/doc/${_cf}"
      )
  endforeach(_cf)
  install(FILES ${_copy_files} DESTINATION "doc")
  list(APPEND _docs ${_copy_files})
else()
  file(GLOB gpin_files RELATIVE "${TOP_SRC_DIR}/doc" "${TOP_SRC_DIR}/doc/*.gp.in")
  set(srcdir "${TOP_SRC_DIR}/doc") # (@srcdir@) will be substituted in *.gp.in
  foreach(_gpin ${gpin_files})
    string(REGEX REPLACE "\\.gp.in$" "" _gp ${_gpin})
    SET_SOURCE_FILES_PROPERTIES("${_gp}.gp" GENERATED)
    configure_file("${TOP_SRC_DIR}/doc/${_gpin}" "${CMAKE_CURRENT_BINARY_DIR}/${_gp}.gp" @ONLY)

    set(_created "${_gp}.txt" "${_gp}.png" "${_gp}.pdf" "${_gp}.eps")
    message(STATUS "About to execute ${GNUPLOT_EXE} ${_gp}.gp")
    add_custom_command(
      OUTPUT ${_created}
      COMMAND ${GNUPLOT_EXE} "${_gp}.gp"
      DEPENDS "${_gp}.gp"
      )
    install(FILES ${_created} DESTINATION "doc")
    list(APPEND _docs ${_created})
  endforeach(_gpin)
endif()

ADD_CUSTOM_TARGET(doc ALL DEPENDS ${_docs})

########### install files ###############




#original Makefile.am contents follow:

#info_TEXINFOS = enblend.texi \
#                enfuse.texi
#
#enblend_TEXINFOS = auxmac.texi auxmac.tex \
#                   fdl.texi \
#                   mask-template-characters.texi \
#                   helpful-programs.texi \
#                   tuning-memory-usage.texi \
#                   understanding-masks.texi \
#                   workflow.texi \
#                   authors.texi \
#                   photographic-workflow.fig \
#                   default.css
#
#enfuse_TEXINFOS = auxmac.texi auxmac.tex \
#                  fdl.texi \
#                  mask-template-characters.texi \
#                  helpful-programs.texi \
#                  tuning-memory-usage.texi \
#                  understanding-masks.texi \
#                  workflow.texi \
#                  authors.texi \
#                  config-edge.gp config.gp \
#                  entropy-cutoff.gp.in \
#                  entropy.gp.in gaussian.gp.in \
#                  laplacian-of-gaussian.gp.in \
#                  sharp-edge.gp.in sharp-edge.data \
#                  smooth-edge.gp.in smooth-edge.data \
#                  local-analysis-window.fig \
#                  photographic-workflow.fig \
#                  focus-stack-decision-tree.fig \
#                  default.css
#
#AM_MAKEINFOFLAGS = @AM_MAKEINFOFLAGS@ \
#                   --css-include=@srcdir@/default.css
#
## created by make(1), user pobably wants to rebuild (often?)
#MOSTLYCLEANFILES =
#
## created by make(1)
#CLEANFILES = entropy.{txt,png,eps,pdf} \
#             entropy-cutoff.{txt,png,eps,pdf} \
#             gaussian.{txt,png,eps,pdf} \
#             laplacian-of-gaussian.{txt,png,eps,pdf} \
#             local-analysis-window.{txt,png,eps,pdf} \
#             photographic-workflow.{txt,png,eps,pdf} \
#             focus-stack-decision-tree.{txt,png,eps,pdf} \
#             sharp-edge.{txt,png,eps,pdf} \
#             smooth-edge.{txt,png,eps,pdf} \
#             *.fig.bak \
#             *.msg
#
## created by configure(1)
#DISTCLEANFILES = texinfo.tex
#
#MAINTAINERCLEANFILES = enblend.info enfuse.info
#
#
## Phony Targets
#
#.PHONY: html-local
#html-local:
#	test -d enblend.html && \
#            for i in $(enblend_TEXINFOS); do \
#              test -f $${i%.png}.png && cp $$i enblend.html; \
#            done; true
#	test -d enfuse.html && \
#            for i in $(enfuse_TEXINFOS); do \
#              test -f $${i%.png}.png && cp $$i enfuse.html; \
#            done; true
#
#.PHONY: clean-local
#clean-local:
#	-rm -rf enblend.xhtml enfuse.xhtml
#
#
## The public id is in "tidy.cfg".
#DTD_SYSID="http://www.w3.org/TR/MathML2/dtd/xhtml-math11-f.dtd"
#
#.PHONY: xhtml
#xhtml: $(HTMLS) html-local
#	for x in $(HTMLS); do \
#	  if test -d $$x; then \
#            test -d $${x/%.html/.xhtml} || mkdir $${x/%.html/.xhtml}; \
#            cp *.png $${x/%.html/.xhtml}; \
#	    for y in $$x/*.html; do \
#	      tidy -config @srcdir@/tidy.cfg -file $${y/%.html/.msg} $$y | \
#	        $(SED) -e '1,9s|""|$(DTD_SYSID)|' \
#                       -e 's|\("[^/"]*\)\.html\([#"]\)|\1.xhtml\2|' > $${y//.html/.xhtml}; \
#            done; \
#	  else \
#	    tidy -config @srcdir@/tidy.cfg -file $${x/%.html/.msg} $$x | \
#	      $(SED) -e '1,9s|""|$(DTD_SYSID)|' > $${x/%.html/.xhtml}; \
#          fi; \
#	done
#
#
## Implicit Rules
#
#.gp.txt:
#	GDFONTPATH=$${GDFONTPATH:-@GDFONTPATH@} $(GNUPLOT) $<
#
#.gp.png:
#	GDFONTPATH=$${GDFONTPATH:-@GDFONTPATH@} $(GNUPLOT) $<
#
#.gp.eps:
#	GDFONTPATH=$${GDFONTPATH:-@GDFONTPATH@} $(GNUPLOT) $<
#
#.gp.pdf:
#	GDFONTPATH=$${GDFONTPATH:-@GDFONTPATH@} $(GNUPLOT) $<
#
#
#.fig.txt:
#	$(SED) -e '1,/---BEGIN-TEXT---/d' \
#               -e '/---END-TEXT---/,$$d' \
#               -e 's/^# \?//' < $< > $@
#
#.fig.png:
#	$(FIG2DEV) -L png $< $@
#
#.fig.eps:
#	$(FIG2DEV) -L eps $< $@
#
#.fig.pdf:
#	$(FIG2DEV) -L pdf $< $@
#
#
## Explicit Rules
#
#enblend.info: photographic-workflow.txt
#
#enblend.dvi: photographic-workflow.eps
#
#enblend.html \
#enblend.pdf: photographic-workflow.png
#
#
#enfuse.info: photographic-workflow.txt \
#             focus-stack-decision-tree.txt \
#             entropy.txt \
#             entropy-cutoff.txt \
#             gaussian.txt \
#             laplacian-of-gaussian.txt \
#             local-analysis-window.txt \
#             sharp-edge.txt \
#             smooth-edge.txt
#
#enfuse.dvi: photographic-workflow.eps \
#            focus-stack-decision-tree.eps \
#            entropy.eps \
#            entropy-cutoff.eps \
#            gaussian.eps \
#            laplacian-of-gaussian.eps \
#            local-analysis-window.eps \
#            sharp-edge.eps \
#            smooth-edge.eps
#
#enfuse.html \
#enfuse.pdf: photographic-workflow.png \
#            focus-stack-decision-tree.png \
#            entropy.png \
#            entropy-cutoff.png \
#            gaussian.png \
#            laplacian-of-gaussian.png \
#            local-analysis-window.png \
#            sharp-edge.png \
#            smooth-edge.png
