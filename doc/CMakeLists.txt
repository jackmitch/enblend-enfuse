# This file is part of enblend.
# Licence details can be found in the file COPYING.
#
# Copyright (c) 2009, Kornel Benko <Kornel.Benko@berlin.de>
#

set(DocumentsDirectory "doc/enblend")
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

macro(copybin _basename)
  set_source_files_properties(${_basename} GENERATED)
  if (${ARGV1} MATCHES "DIR")
    set(cp_command "copy_directory")
    install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${_basename}" DESTINATION ${DocumentsDirectory})
  else()
    set(cp_command "copy")
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${_basename}" DESTINATION ${DocumentsDirectory})
  endif()
  add_custom_command(
    OUTPUT ${_basename}
    COMMAND ${CMAKE_COMMAND} -E ${cp_command} "${TOP_SRC_DIR}/doc/${_basename}" ${_basename}
    DEPENDS "${TOP_SRC_DIR}/doc/${_basename}"
    )
  list(APPEND _docs ${_basename})
  list(APPEND _depend_gp ${_basename})
endmacro(copybin)

project(doc1)

# search for gnuplot and
# generate images (.txt, .pdf .eps and .png)
find_program(GNUPLOT_EXE "gnuplot")
set(_depend_gp)
if (NOT DOC OR ${GNUPLOT_EXE} MATCHES "-NOTFOUND")
  # there is no gnuplot command, so simply copy the needed files
  file(GLOB _copy_files RELATIVE "${TOP_SRC_DIR}/doc"
    "${TOP_SRC_DIR}/doc/*.txt"
    "${TOP_SRC_DIR}/doc/*.eps"
    "${TOP_SRC_DIR}/doc/*.png"
    "${TOP_SRC_DIR}/doc/*.pdf")
  foreach(_cf ${_copy_files})
    if(NOT ${_cf} MATCHES "en(blend|fuse)") # This is done below
      copybin("${_cf}")
    endif()
  endforeach(_cf)
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
    foreach(_cr ${_created})
      install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${_cr}" DESTINATION ${DocumentsDirectory})
    endforeach()
    list(APPEND _docs ${_created})
    list(APPEND _depend_gp ${_created})
  endforeach(_gpin)
endif()

project(doc2)
# Search for fig2dev and convert some images
find_program(FIG2DEV_EXE fig2dev)
if (${FIG2DEV_EXE} MATCHES "-NOTFOUND")
  # here the creation of html and pdf will be disabled, because .gif-prerequisites are not found
  # Automake uses here sed and the like
  set(DOC3:BOOL "FALSE")
else()
  foreach(_img "photographic-workflow" "focus-stack-decision-tree" "local-analysis-window")
    set(_created "${CMAKE_CURRENT_BINARY_DIR}/${_img}")
    foreach(_ext gif pdf)
      add_custom_command(
	OUTPUT "${_created}.${_ext}"
	COMMAND ${FIG2DEV_EXE} -L ${_ext} "${TOP_SRC_DIR}/doc/${_img}.fig" ">" "${_created}.${_ext}"
	DEPENDS "${TOP_SRC_DIR}/doc/${_img}.fig"
	)
      list(APPEND _docs ${_created}.${_ext})
      list(APPEND _depend_gp ${_created}.${_ext})
    endforeach(_ext)
  endforeach(_img)
  set(DOC3:BOOL ${DOC})
endif()

project(doc3)
# Search for makeinfo command and
# generate info and pdf for enblend + enfuse
find_program(MAKEINFO_EXE "makeinfo")
find_program(TEXI2DVI_EXE "texi2dvi")
foreach (_cmd "enblend" "enfuse")
  set(_created "${CMAKE_CURRENT_BINARY_DIR}/${_cmd}")
  if (NOT DOC3:BOOL OR ${MAKEINFO_EXE} MATCHES "-NOTFOUND" OR ${TEXI2DVI_EXE} MATCHES "-NOTFOUND")
    # there is no makeinfo or no texi2dvi command,
    # this implies insufficient latex installation.
    # Therefore simply copy the data from source tree.
    foreach(_ext "info")
      copybin("${_cmd}.${_ext}")
    endforeach()
    #
    # Disabled, because not part of repo.
    # They would have only been created with "make pdf" or "make html" in source doc directory
    # under automake
    if (EXISTS "${TOP_SRC_DIR}/doc/${_cmd}.pdf" AND NOT IS_DIRECTORY "${TOP_SRC_DIR}/doc/${_cmd}.pdf")
      copybin("${_cmd}.pdf")
    else()
      message(STATUS "Sorry, but there is no file ${TOP_SRC_DIR}/doc/${_cmd}.pdf to copy. Leave it out")
    endif()
    if (EXISTS "${TOP_SRC_DIR}/doc/${_cmd}.html" AND IS_DIRECTORY "${TOP_SRC_DIR}/doc/${_cmd}.html")
      copybin("${_cmd}.html" "DIR")
    else()
      message(STATUS "Sorry, but there is no directory ${TOP_SRC_DIR}/doc/${_cmd}.html to copy. Leave it out")
    endif()
  else()
    # here fig2dev, makeinfo and texi2dvi exist, option DOC set to ON
    message(STATUS "create ${_created}.info")
    add_custom_command(
      OUTPUT "${_created}.info"
      COMMAND ${MAKEINFO_EXE} "--css-include=${TOP_SRC_DIR}/doc/default.css"
      "-I" "${TOP_SRC_DIR}/doc"
      -o "${_created}.info" "${TOP_SRC_DIR}/doc/${_cmd}.texi"
      DEPENDS "${TOP_SRC_DIR}/doc/${_cmd}.texi" ${_${_cmd}_TEXINFOS} ${_depend_gp}
      )
    install(FILES "${_created}.info" DESTINATION ${DocumentsDirectory})
    add_custom_command(
      OUTPUT "${_created}.pdf"
      COMMAND ${TEXI2DVI_EXE} "--pdf" "--batch" "${TOP_SRC_DIR}/doc/${_cmd}.texi"
      "${TOP_SRC_DIR}/doc/${_cmd}.texi"
      DEPENDS "${TOP_SRC_DIR}/doc/${_cmd}.texi" ${_${_cmd}_TEXINFOS} ${_depend_gp}
      )
    install(FILES "${_created}.pdf" DESTINATION ${DocumentsDirectory})
    message(STATUS "create ${_created}.html")
    add_custom_command(
      OUTPUT "${_created}.html/index.html"
      COMMAND ${MAKEINFO_EXE} "--html" "--css-include=${TOP_SRC_DIR}/doc/default.css"
      "-I" "${TOP_SRC_DIR}/doc"
      -o "${_created}.html" "${TOP_SRC_DIR}/doc/${_cmd}.texi"
      DEPENDS "${TOP_SRC_DIR}/doc/${_cmd}.texi" ${_${_cmd}_TEXINFOS} ${_depend_gp}
      )
    install(DIRECTORY "${_created}.html" DESTINATION ${DocumentsDirectory})
  endif()
  list(APPEND _docs "${_created}.info" "${_created}.pdf" "${_created}.html")
endforeach(_cmd)

add_custom_target(doc3 ALL DEPENDS ${_docs})

# set make order
add_dependencies(doc3 doc2 doc1 enblendenfuse)


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
