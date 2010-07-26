# This file is part of enblend.
# Licence details can be found in the file COPYING.
#
# Copyright (c) 2009-2010, Kornel Benko <Kornel.Benko@berlin.de>
#

include(GetPrerequisites)	# get function gp_append_unique()
IF(WIN32)
set(DocumentsDirectory "doc")
ELSE(WIN32)
set(DocumentsDirectory "doc/enblend")
ENDIF(WIN32)

set(_docs)
set(_enblend_Bases auxmac.texi auxmac.tex fdl.texi mask-template-characters.texi helpful-programs.texi
    tuning-memory-usage.texi understanding-masks.texi workflow.texi authors.texi color-profiles.texi
    photographic-workflow.fig default.css)

set(_enfuse_Bases auxmac.texi auxmac.tex fdl.texi mask-template-characters.texi helpful-programs.texi
    tuning-memory-usage.texi understanding-masks.texi workflow.texi authors.texi color-profiles.texi
    config-edge.gp config.gp entropy-cutoff.gp.in entropy.gp.in gaussian.gp.in
    laplacian-of-gaussian.gp.in sharp-edge.gp.in sharp-edge.data
    smooth-edge.gp.in smooth-edge.data local-analysis-window.fig
    photographic-workflow.fig focus-stack-decision-tree.fig default.css)

set(_enblend_TEXINFOS)
foreach(_e ${_enblend_Bases})
  list(APPEND _enblend_TEXINFOS "${TOP_SRC_DIR}/doc/${_e}")
endforeach(_e)

set(_enfuse_TEXINFOS)
foreach(_e ${_enfuse_Bases})
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
    COMMENT "Copy ${TOP_SRC_DIR}/doc/${_basename} ${_basename}"
    )
  list(APPEND _docs ${_basename})
  list(APPEND _depend_gp ${_basename})
endmacro(copybin)

macro(create_selection _list _ext _destlist)
  set(${_destlist})		# empty destination list first
  foreach(_src ${${_list}})
    if(${_src} MATCHES "^(.+)\\.${_ext}$") #select from the list those with expected extension
      list(APPEND ${_destlist} ${CMAKE_MATCH_1}) # Add the base name without extension only
    endif()
  endforeach(_src)
endmacro(create_selection)

set(_enblendenfuse_fig_images)
set(_enblendenfuse_gp_images)

foreach(_cmd "enblend" "enfuse")
  create_selection(_${_cmd}_Bases "fig" _${_cmd}_fig_images)
  create_selection(_${_cmd}_Bases "gp\\\\.in" _${_cmd}_gp_images)
  set(_${_cmd}_images ${_${_cmd}_fig_images} ${_${_cmd}_gp_images})
  foreach(_i ${_${_cmd}_fig_images})
    gp_append_unique(_enblendenfuse_fig_images ${_i})
  endforeach(_i)
  foreach(_i ${_${_cmd}_gp_images})
    gp_append_unique(_enblendenfuse_gp_images ${_i})
  endforeach(_i)
endforeach(_cmd)

project(doc0)

set(_helpfiles)

#create config.h.texi
add_custom_command(
  OUTPUT  "${CMAKE_CURRENT_BINARY_DIR}/config-h.texi"
  COMMAND ${PERL_EXECUTABLE} "${TOP_SRC_DIR}/doc/define2set.pl"
  "${CMAKE_BINARY_DIR}/config.h" ">" "${CMAKE_CURRENT_BINARY_DIR}/config-h.texi"
  DEPENDS "${CMAKE_BINARY_DIR}/config.h"
  COMMENT "Creating config-h.texi"
  )
SET_SOURCE_FILES_PROPERTIES("${CMAKE_CURRENT_BINARY_DIR}/config-h.texi" GENERATED)
list(APPEND _helpfiles "${CMAKE_CURRENT_BINARY_DIR}/config-h.texi")

# create versenblend.texi and versenfuse.texi
foreach(_cmd enblend enfuse)
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/vers${_cmd}.texi"
    COMMAND ${PERL_EXECUTABLE} "${TOP_SRC_DIR}/doc/CreateVersTexi.pl" "${TOP_SRC_DIR}/doc/${_cmd}.texi" "${_version_lines}" ${A4} > "${CMAKE_CURRENT_BINARY_DIR}/vers${_cmd}.texi"
    DEPENDS "${TOP_SRC_DIR}/doc/${_cmd}.texi"
    COMMENT "perl: Creating vers${_cmd}.texi"
    )
  SET_SOURCE_FILES_PROPERTIES("${CMAKE_CURRENT_BINARY_DIR}/vers${_cmd}.texi" GENERATED)
  list(APPEND _helpfiles "${CMAKE_CURRENT_BINARY_DIR}/vers${_cmd}.texi")
endforeach(_cmd)

# create varsenblend.texi and varsenfuse.texi
foreach(_cmd enblend enfuse)
  set(_depends "${TOP_SRC_DIR}/src/${_cmd}.cc"
    "${TOP_SRC_DIR}/src/bounds.h" "${TOP_SRC_DIR}/src/common.h"
    "${TOP_SRC_DIR}/src/global.h")
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/vars${_cmd}.texi"
    COMMAND ${PERL_EXECUTABLE} "${TOP_SRC_DIR}/doc/docstrings" ${_depends} ">" vars${_cmd}.texi
    DEPENDS ${_depends}
    COMMENT "perl: creating vars${_cmd}.texi"
    )
  SET_SOURCE_FILES_PROPERTIES("${CMAKE_CURRENT_BINARY_DIR}/vars${_cmd}.texi" GENERATED)
  list(APPEND _helpfiles "${CMAKE_CURRENT_BINARY_DIR}/vars${_cmd}.texi")
endforeach(_cmd)

foreach(_fig ${_enblendenfuse_fig_images})
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${_fig}.txt"
    COMMAND ${PERL_EXECUTABLE} "${TOP_SRC_DIR}/doc/fig2txt.pl" "${TOP_SRC_DIR}/doc/${_fig}.fig" ">" "${CMAKE_CURRENT_BINARY_DIR}/${_fig}.txt"
    DEPENDS "${TOP_SRC_DIR}/doc/${_fig}.fig"
    COMMENT "perl: creating ${CMAKE_CURRENT_BINARY_DIR}/${_fig}.txt"
    )
  SET_SOURCE_FILES_PROPERTIES("${CMAKE_CURRENT_BINARY_DIR}/${_fig}.txt" GENERATED)
  list(APPEND _helpfiles "${CMAKE_CURRENT_BINARY_DIR}/${_fig}.txt")
endforeach(_fig)
add_custom_target(doc0 ALL DEPENDS ${_helpfiles})

project(doc1)

# search for gnuplot and
# generate images (.txt, .pdf .eps and .png)
find_package(Gnuplot)
set(_depend_gp)
if (NOT DOC OR NOT GNUPLOT_FOUND)
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
  #
  set(srcdir "${TOP_SRC_DIR}/doc") # (@srcdir@) will be substituted in *.gp.in
  foreach(_gp ${_enblendenfuse_gp_images})
    configure_file("${TOP_SRC_DIR}/doc/${_gp}.gp.in" "${CMAKE_CURRENT_BINARY_DIR}/${_gp}.gp" @ONLY)

    set(_created "${CMAKE_CURRENT_BINARY_DIR}/${_gp}.txt" "${CMAKE_CURRENT_BINARY_DIR}/${_gp}.png" "${CMAKE_CURRENT_BINARY_DIR}/${_gp}.pdf" "${CMAKE_CURRENT_BINARY_DIR}/${_gp}.eps")
    #message(STATUS "About to execute ${GNUPLOT_EXECUTABLE} ${_gp}.gp")
    add_custom_command(
      OUTPUT ${_created}
      COMMAND ${GNUPLOT_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/${_gp}.gp"
      DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${_gp}.gp"
      COMMENT "gnuplot: Creating ${_created} from ${CMAKE_CURRENT_BINARY_DIR}/${_gp}.gp"
      )
    SET_SOURCE_FILES_PROPERTIES(${_created} GENERATED)
    install(FILES ${_created} DESTINATION ${DocumentsDirectory})
    #message(STATUS "Generated files ${_created}")
    list(APPEND _docs ${_created})
    list(APPEND _depend_gp ${_created})
  endforeach(_gp)
endif()
add_custom_target(doc1 ALL DEPENDS ${_depend_gp})
add_dependencies(doc1 doc0)
project(doc2)

# Search for fig2dev and convert some images
find_program(FIG2DEV_EXE fig2dev)
if (${FIG2DEV_EXE} MATCHES "-NOTFOUND")
  # here the creation of html and pdf will be disabled, because .gif-prerequisites are not found
  # Automake uses here sed and the like
  set(DOC3:BOOL "FALSE")
else()
  #message(STATUS "_enblendenfuse_fig_images =  ${_enblendenfuse_fig_images}")
  foreach(_img ${_enblendenfuse_fig_images})
    set(_created "${CMAKE_CURRENT_BINARY_DIR}/${_img}")
    foreach(_ext gif pdf png)
      #message(STATUS "will create ${_created}.${_ext}")
      add_custom_command(
	OUTPUT "${_created}.${_ext}"
	COMMAND ${FIG2DEV_EXE} -L ${_ext} "${TOP_SRC_DIR}/doc/${_img}.fig" ">" "${_created}.${_ext}"
	DEPENDS "${TOP_SRC_DIR}/doc/${_img}.fig"
	COMMENT "fig2dev: Creating ${_created}.${_ext} from ${TOP_SRC_DIR}/doc/${_img}.fig"
	)
      list(APPEND _docs ${_created}.${_ext})
      list(APPEND _depend_gp ${_created}.${_ext})
    endforeach(_ext)
    set_source_files_properties(${_created}.png GENERATED)
    if(NOSPLIT)
      install(FILES ${_created}.png DESTINATION ${DocumentsDirectory})
    endif(NOSPLIT)
  endforeach(_img)
  set(DOC3:BOOL ${DOC})
endif()
add_custom_target(doc2 ALL DEPENDS ${_depend_gp})
add_dependencies(doc2 doc1)

project(doc3)
# Search for makeinfo command and
# generate info and pdf for enblend + enfuse
find_program(MAKEINFO_EXE "makeinfo")
if (DOC3:BOOL)
  if ( NOT ${MAKEINFO_EXE} MATCHES "-NOTFOUND")
    # Check for sufficient version, should be >= 4.11
    EXECUTE_PROCESS(COMMAND ${MAKEINFO_EXE} "--version" OUTPUT_VARIABLE _makeinfo_version OUTPUT_STRIP_TRAILING_WHITESPACE)
    foreach(_v_l ${_makeinfo_version})
      #message(STATUS "_v_l = ${_v_l}")
      if(_v_l MATCHES "\\(GNU texinfo\\) *\([0-9]+\)\\.\([0-9]+\)")
	set(MAKEINFO_MAJOR ${CMAKE_MATCH_1})
	set(MAKEINFO_MINOR1 "00000${CMAKE_MATCH_2}")
	if(MAKEINFO_MINOR1 MATCHES "^[0]+\([0-9][0-9][0-9][0-9]\)\$")
	  set(MAKEINFO_MINOR ${CMAKE_MATCH_1})
	endif()
	set(MAKEINFO_VERSION "${MAKEINFO_MAJOR}.${MAKEINFO_MINOR}")
      endif()
    endforeach()
    string(COMPARE GREATER "${MAKEINFO_VERSION}" "4.0010" DOC3:BOOL)
    #message(STATUS "MAKEINFO_VERSION = ${MAKEINFO_VERSION}")
    if (NOT DOC3:BOOL)
      message(STATUS "Sorry, but you need makeinfo version >= 4.11")
    endif()
  endif()
endif()

find_program(TEXI2DVI_EXE "texi2dvi")
if(${EXTRATEXI2DVIFLAGS} MATCHES "^\\@")
  set(TEXI2DVIFLAGS "--command=${EXTRATEXI2DVIFLAGS}")
  message(STATUS "EXTRATEXI2DVIFLAGS = ${EXTRATEXI2DVIFLAGS}")
  message(STATUS "Setting TEXI2DVIFLAGS to ${TEXI2DVIFLAGS}")
else()
  set(TEXI2DVIFLAGS)
  message(STATUS "Resetting TEXI2DVIFLAGS")
endif()

set(_auxdepends)
foreach(_f auxmac.texi auxmac.tex workflow.texi filespec.texi mask-template-characters.texi
    understanding-masks.texi tuning-memory-usage.texi helpful-programs.texi bug-reports.texi
    authors.texi fdl.texi color-profiles.texi)
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${_f}"
    COMMAND ${CMAKE_COMMAND} -E copy "${TOP_SRC_DIR}/doc/${_f}" "${CMAKE_CURRENT_BINARY_DIR}/${_f}"
    DEPENDS "${TOP_SRC_DIR}/doc/${_f}"
    )
  list(APPEND _auxdepends "${CMAKE_CURRENT_BINARY_DIR}/${_f}")
endforeach()

foreach (_cmd "enblend" "enfuse")
  set(_created "${CMAKE_CURRENT_BINARY_DIR}/${_cmd}")
  if (NOT DOC3:BOOL OR ${MAKEINFO_EXE} MATCHES "-NOTFOUND" OR ${TEXI2DVI_EXE} MATCHES "-NOTFOUND")
    # there is no makeinfo or no texi2dvi command,
    # this implies insufficient latex installation.
    # Therefore simply copy the data from source tree.
    foreach(_ext "info")
      copybin("${_cmd}.${_ext}")
      list(APPEND _docs "${_created}.${_ext}")
    endforeach()
    #
    # Disabled, because not part of repo.
    # They would have only been created with "make pdf" or "make html" in source doc directory
    # under automake
    if (EXISTS "${TOP_SRC_DIR}/doc/${_cmd}.pdf" AND NOT IS_DIRECTORY "${TOP_SRC_DIR}/doc/${_cmd}.pdf")
      copybin("${_cmd}.pdf")
      list(APPEND _docs "${_created}.pdf")
    else()
      message(STATUS "Sorry, but there is no file ${TOP_SRC_DIR}/doc/${_cmd}.pdf to copy. Leave it out")
    endif()
    if (EXISTS "${TOP_SRC_DIR}/doc/${_cmd}.html" AND IS_DIRECTORY "${TOP_SRC_DIR}/doc/${_cmd}.html")
      copybin("${_cmd}.html" "DIR")
      list(APPEND _docs "${_created}.html")
    else()
      message(STATUS "Sorry, but there is no directory ${TOP_SRC_DIR}/doc/${_cmd}.html to copy. Leave it out")
    endif()
  else()
    # here fig2dev, makeinfo and texi2dvi exist, option DOC set to ON
    #message(STATUS "create ${_created}.info")
    add_custom_command(
      OUTPUT "${_created}.info"
      COMMAND ${MAKEINFO_EXE} "--css-include=${TOP_SRC_DIR}/doc/default.css"
      "-I" "${TOP_SRC_DIR}/doc"
      -o "${_created}.info" "${TOP_SRC_DIR}/doc/${_cmd}.texi"
      DEPENDS "${TOP_SRC_DIR}/doc/${_cmd}.texi"
      ${_${_cmd}_TEXINFOS}
      ${_depend_gp}
      "${CMAKE_CURRENT_BINARY_DIR}/config-h.texi"
      COMMENT "makeinfo: Creating ${_created}.info from ${TOP_SRC_DIR}/doc/${_cmd}.texi"
      )
    install(FILES "${_created}.info" DESTINATION ${DocumentsDirectory})
    add_custom_command(
      OUTPUT "${_created}.texi"
      COMMAND ${CMAKE_COMMAND} -E copy "${TOP_SRC_DIR}/doc/${_cmd}.texi" "${_created}.texi"
      DEPENDS "${TOP_SRC_DIR}/doc/${_cmd}.texi"
      )
    add_custom_command(
      OUTPUT "${_created}.pdf"
      COMMAND ${TEXI2DVI_EXE} ${TEXI2DVIFLAGS} "--pdf" "--batch" "${_created}.texi"
      DEPENDS "${_created}.texi"
      ${_${_cmd}_TEXINFOS} ${_depend_gp}
      "${CMAKE_CURRENT_BINARY_DIR}/config-h.texi" ${_auxdepends}
      COMMENT "texi2dvi: Creating ${_created}.pdf from ${TOP_SRC_DIR}/doc/${_cmd}.texi, cmdline="${TEXI2DVIFLAGS}
      )
    install(FILES "${_created}.pdf" DESTINATION ${DocumentsDirectory})
    #message(STATUS "create ${_created}.html")
    if(NOSPLIT)
      set(_created_output "${_created}_file.html")
      set(_output_param "${_created}_file.html")
      set(_html_param "--no-split")
    else()
      set(_created_output "${_created}.html/index.html")
      set(_output_param "${_created}.html")
      set(_html_param)
    endif(NOSPLIT)
    add_custom_command(
      OUTPUT "${_created_output}"
      COMMAND ${MAKEINFO_EXE} "--html" ${_html_param} "--css-include=${TOP_SRC_DIR}/doc/default.css"
      "-I" "${TOP_SRC_DIR}/doc"
      -o "${_output_param}" "${TOP_SRC_DIR}/doc/${_cmd}.texi"
      DEPENDS "${TOP_SRC_DIR}/doc/${_cmd}.texi" ${_${_cmd}_TEXINFOS} ${_depend_gp}
      COMMENT "makeinfo: Creating ${_created}.html from ${TOP_SRC_DIR}/doc/${_cmd}.texi"
      )
    if(NOSPLIT)
      install(FILES "${_created}_file.html" DESTINATION ${DocumentsDirectory} RENAME "${_cmd}.html")
    else(NOSPLIT)
      foreach(_img ${_${_cmd}_images})
	set(_needed_image "${_created}.html/${_img}.png")
	set(_parent_image "${CMAKE_CURRENT_BINARY_DIR}/${_img}.png")
	#message(STATUS "_needed_image = ${_needed_image}")
	#message(STATUS "_parent_image = ${_parent_image}")
	add_custom_command(
	  OUTPUT "${_needed_image}"
	  COMMAND ${CMAKE_COMMAND} -E copy "${_parent_image}" "${_needed_image}"
	  DEPENDS "${_created}.html/index.html" "${_parent_image}"
	  COMMENT "Copying ${_parent_image} to ${_needed_image}"
	  )
	list(APPEND _docs "${_needed_image}")
      endforeach(_img)
      install(DIRECTORY "${_created}.html" DESTINATION ${DocumentsDirectory})
    endif(NOSPLIT)
    list(APPEND _docs "${_created}.info" "${_created}.pdf" "${_created_output}")
  endif()
endforeach(_cmd)

add_custom_target(doc3 ALL DEPENDS ${_docs})

if(DOC3:BOOL)
  # set make order
  add_dependencies(doc3 doc2 doc1 doc0 enblendenfuse)
endif()


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
