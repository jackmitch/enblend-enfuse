
project(doc)

IF(WIN32)
set(DocumentsDirectory "doc")
ELSE(WIN32)
set(DocumentsDirectory "share/doc/enblend")
ENDIF(WIN32)

# some used flags
set(UPDATED_ON_FLAGS "--database=source-code-manager --override=none:\\datebyunknown --override=fs:\\datebyfs --override=scm:\\datebyscm")

set(GNUPLOT_FLAGS "--default-settings -e DATA_DIR='${CMAKE_CURRENT_SOURCE_DIR}'")
set(GNUPLOT_PNG_TERMINAL_OPTIONS "enhanced transparent truecolor size 720,540")
set(GNUPLOT_EPS_TERMINAL_OPTIONS "color size 4.5in,3.375in")

set(DOTPP_FLAGS "--fatal-warnings --prefix-builtins --synclines")
set(DOTPP_PRELOAD "${CMAKE_CURRENT_SOURCE_DIR}/uml-dot.m4")
set(DOT_FLAGS "")
set(DOT_EPS_FLAGS "-Gsize=4.75,8.5")
set(DOT_PNG_FLAGS "-Gbgcolor=transparent -Gresolution=${PRINTRESOLUTION}")

set(SVGCONVERT_FLAGS "--dpi-x=96 --dpi-y=96")
set(SVGCONVERT_EPS_FLAGS "${SVGCONVERT_FLAGS} --format=eps")
set(SVGCONVERT_PNG_FLAGS "${SVGCONVERT_FLAGS} --format=png")

set(LATEX_FLAGS "-file-line-error -halt-on-error -interaction nonstopmode")
set(PDFLATEX_FLAGS "-file-line-error -halt-on-error -interaction nonstopmode")
set(MAKEINDEX_FLAGS "")
set(DVIPS_FLAGS "-q -t ${PAPERSIZE} -D ${PRINTRESOLUTION} -Z")

# Some other "interesting" flags for Hevea:
#     -dv        add borders to block-level elements
#     -mathml    enable MathML output
set(LATEX2HTML_FLAGS "-O -I . -I ${CMAKE_CURRENT_SOURCE_DIR} -fix")

# format arguments for later usage
separate_arguments(UPDATED_ON_FLAGS)
separate_arguments(GNUPLOT_FLAGS)
separate_arguments(DOTPP_FLAGS)
separate_arguments(DOT_FLAGS)
separate_arguments(DOT_EPS_FLAGS)
separate_arguments(DOT_PNG_FLAGS)
separate_arguments(SVGCONVERT_EPS_FLAGS)
separate_arguments(SVGCONVERT_PNG_FLAGS)
separate_arguments(LATEX_FLAGS)
separate_arguments(PDFLATEX_FLAGS)
separate_arguments(MAKEINDEX_FLAGS)
separate_arguments(DVIPS_FLAGS)
separate_arguments(LATEX2HTML_FLAGS)

# build some dynamic files needed for doc
add_custom_command(
    OUTPUT  "${CMAKE_BINARY_DIR}/doc/config-h.tex"
    COMMAND ${PERL_EXECUTABLE} "${TOP_SRC_DIR}/doc/config2tex" "${CMAKE_BINARY_DIR}/config.h" ">" "${CMAKE_BINARY_DIR}/doc/config-h.tex"
    DEPENDS "${CMAKE_BINARY_DIR}/config.h"
)

set(_dynamic_preamble "\\PassOptionsToPackage{${PAPERSIZE}paper}{report}\n\\PassOptionsToPackage{${PAPERSIZE}paper}{refrep}\n\\newif\\ifreferencemanual\n\\newif\\ifhyperref\n${DYNAMIC_TEX_PREAMBLE}\n")
if(INSTALL_PDF_DOC AND NOT ${_dynamic_preamble} MATCHES "\\hyperreftrue")
  string(REPLACE "\\newif\\ifhyperref\n" "\\newif\\ifhyperref\n\\hyperreftrue\n" _dynamic_preamble "${_dynamic_preamble}")
endif()
file(WRITE "${CMAKE_BINARY_DIR}/doc/dynamic-preamble.tex" "${_dynamic_preamble}")
set_source_files_properties("${CMAKE_BINARY_DIR}/doc/dynamic-preamble.tex" PROPERTIES GENERATED TRUE)

set(DOCSTRING_SOURCE_FILES
    "${TOP_SRC_DIR}/src/bounds.h"
    "${TOP_SRC_DIR}/src/common.h"
    "${TOP_SRC_DIR}/src/filespec.h"
    "${TOP_SRC_DIR}/src/filespec.cc"
    "${TOP_SRC_DIR}/src/global.h"
    "${TOP_SRC_DIR}/src/mga.h"
    "${TOP_SRC_DIR}/src/opencl.cc"
)
add_custom_command(
    OUTPUT  "${CMAKE_BINARY_DIR}/doc/enblend-variables.tex"
    COMMAND ${PERL_EXECUTABLE} "${TOP_SRC_DIR}/doc/docstrings" "scan-files" ${DOCSTRING_SOURCE_FILES} "${TOP_SRC_DIR}/src/enblend.cc" ">" "${CMAKE_BINARY_DIR}/doc/enblend-variables.tex"
    DEPENDS ${DOCSTRING_SOURCE_FILES} "${TOP_SRC_DIR}/src/enblend.cc"
    COMMENT "Creating enblend-variables.tex"
)

add_custom_command(
    OUTPUT  "${CMAKE_BINARY_DIR}/doc/enfuse-variables.tex"
    COMMAND ${PERL_EXECUTABLE} "${TOP_SRC_DIR}/doc/docstrings" "scan-files" ${DOCSTRING_SOURCE_FILES} "${TOP_SRC_DIR}/src/enfuse.cc" ">" "${CMAKE_BINARY_DIR}/doc/enfuse-variables.tex"
    DEPENDS ${DOCSTRING_SOURCE_FILES} "${TOP_SRC_DIR}/src/enfuse.cc"
    COMMENT "Creating enfuse-variables.tex"
)

if(WIN32)
  SET(INSTALL_SHARE "share")
else()
  SET(INSTALL_SHARE "${CMAKE_INSTALL_PREFIX}/$&%_/share")
endif()
# replace invalid characters for latex
string(REGEX REPLACE "([$&%_^{}])" "\\\\\\1" INSTALL_SHARE "${INSTALL_SHARE}")
configure_file(version.cmake.in "${CMAKE_BINARY_DIR}/doc/version.in" @ONLY)

file(GLOB ENBLEND_SOURCE_FILES common-*.tex enblend-*.tex enblend.tex static-preamble.tex lead-in.tex lead-out.tex index-def.tex index-use.tex
              macros.tex exemplar.hva htmlstyle.tex url-def.tex
)
add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/doc/enblend-version.tex" "${CMAKE_BINARY_DIR}/doc/enblend_date"
    COMMAND ${PERL_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/updated-on ${UPDATED_ON_FLAGS} ${ENBLEND_SOURCE_FILES} > ${CMAKE_BINARY_DIR}/doc/enblend_date
    COMMAND ${CMAKE_COMMAND} -DINPUT_PATH=${CMAKE_BINARY_DIR}/doc -DFILENAME=enblend-version.tex -DDATEFILE=enblend_date -P "${CMAKE_CURRENT_SOURCE_DIR}/updateDate.cmake"
    DEPENDS "${CMAKE_BINARY_DIR}/doc/version.in" ${ENBLEND_SOURCE_FILES}
    VERBATIM
)
set(ENBLEND_GENERATED_FILES
    config-h.tex
    dynamic-preamble.tex
    enblend-variables.tex
    enblend-version.tex
)

file(GLOB ENFUSE_SOURCE_FILES common-*.tex enfuse-*.tex enfuse.tex static-preamble.tex lead-in.tex lead-out.tex index-def.tex index-use.tex
              macros.tex exemplar.hva htmlstyle.tex url-def.tex
)
add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/doc/enfuse-version.tex" "${CMAKE_BINARY_DIR}/doc/enfuse_date"
    COMMAND ${PERL_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/updated-on ${UPDATED_ON_FLAGS} ${ENFUSE_SOURCE_FILES} > ${CMAKE_BINARY_DIR}/doc/enfuse_date
    COMMAND ${CMAKE_COMMAND} -DINPUT_PATH="${CMAKE_BINARY_DIR}/doc" -DFILENAME=enfuse-version.tex -DDATEFILE=enfuse_date -P "${CMAKE_CURRENT_SOURCE_DIR}/updateDate.cmake"
    DEPENDS "${CMAKE_BINARY_DIR}/doc/version.in" ${ENFUSE_SOURCE_FILES}
)
set(ENFUSE_GENERATED_FILES
    config-h.tex
    dynamic-preamble.tex
    enfuse-variables.tex
    enfuse-version.tex
)

add_custom_target(generated_files_for_doc DEPENDS ${ENBLEND_GENERATED_FILES} ${ENFUSE_GENERATED_FILES})

# prepare picture for ps/pdf version
if(WIN32)
    # MiKTeX epstopdf does not understand PostscriptLevel 3 file, so fall back to normal eps
    add_custom_command(
        OUTPUT "${CMAKE_BINARY_DIR}/doc/seam-line-visualization.eps"
        COMMAND ${IMAGEMAGICK_CONVERT_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/seam-line-visualization.tif "eps:${CMAKE_BINARY_DIR}/doc/seam-line-visualization.eps"
        DEPENDS seam-line-visualization.tif
    )
else()
    add_custom_command(
        OUTPUT "${CMAKE_BINARY_DIR}/doc/seam-line-visualization.eps"
        COMMAND ${IMAGEMAGICK_CONVERT_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/seam-line-visualization.tif "eps3:${CMAKE_BINARY_DIR}/doc/seam-line-visualization.eps"
        DEPENDS seam-line-visualization.tif
    )
endif()
set(ENBLEND_PSTEXFILES "${CMAKE_BINARY_DIR}/doc/seam-line-visualization.eps")

macro(build_pstex _list _depends)
    foreach(_file ${_list})
        get_filename_component(_name ${_file} NAME_WE)
        get_filename_component(_ext ${_file} EXT)
        if(_ext STREQUAL ".svg")
            add_custom_command(
                OUTPUT "${CMAKE_BINARY_DIR}/doc/${_name}.eps"
                COMMAND ${RSVGCONVERT_EXECUTABLE} ${SVGCONVERT_EPS_FLAGS} --output=${_name}.eps ${CMAKE_CURRENT_SOURCE_DIR}/${_name}.svg
                DEPENDS ${_name}.svg
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc
                VERBATIM
            )
            list(APPEND ${_depends} "${CMAKE_BINARY_DIR}/doc/${_name}.eps")
        elseif(_ext STREQUAL ".dot")
            add_custom_command(
                OUTPUT "${CMAKE_BINARY_DIR}/doc/${_name}.eps"
                COMMAND ${M4_EXECUTABLE} ${DOTPP_FLAGS} "--define=dot_output_type=eps" ${DOTPP_PRELOAD} ${CMAKE_CURRENT_SOURCE_DIR}/${_name}.dot | ${DOT_EXECUTABLE} ${DOT_FLAGS} -Teps ${DOT_EPS_FLAGS} -o${_name}.eps
                DEPENDS ${_name}.dot ${CMAKE_CURRENT_SOURCE_DIR}/uml-dot.m4
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc
                VERBATIM
            )
            list(APPEND ${_depends} "${CMAKE_BINARY_DIR}/doc/${_name}.eps")
        elseif(_ext STREQUAL ".gp")
            add_custom_command(
                OUTPUT "${CMAKE_BINARY_DIR}/doc/${_name}.eps"
                COMMAND ${GNUPLOT_EXECUTABLE} ${GNUPLOT_FLAGS} -e "set output \"${_name}.eps\"; set terminal eps ${GNUPLOT_EPS_TERMINAL_OPTIONS}" ${CMAKE_CURRENT_SOURCE_DIR}/${_name}.gp
                DEPENDS ${_name}.gp colors.gp
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc
                VERBATIM
            )
            list(APPEND ${_depends} "${CMAKE_BINARY_DIR}/doc/${_name}.eps")
        else()
            # this should never happens
            message(FATAL_ERROR "Unknown extension")
        endif()
    endforeach()
endmacro()

set(ENBLEND_IMAGE_INPUTFILES
    external-mask-workflow.dot
    photographic-workflow.dot
    internal-enblend-flow.dot
    internal-optimizer-chain.dot
    log-transform.gp
)
build_pstex("${ENBLEND_IMAGE_INPUTFILES}" ENBLEND_PSTEXFILES)

set(ENFUSE_PSTEXFILES)
set(ENFUSE_IMAGE_INPUTFILES
    external-mask-workflow.dot
    focus-stack-decision-tree.dot
    local-analysis-window.svg
    photographic-workflow.dot
    entropy-cutoff.gp
    entropy.gp
    exposure-cutoff.gp
    exposure-weights.gp
    fullsine.gp
    gaussian.gp
    halfsine.gp
    laplacian-of-gaussian.gp
    log-transform.gp
    lorentzian.gp
    power.gp
    sharp-edge.gp
    smooth-edge.gp
)
build_pstex("${ENFUSE_IMAGE_INPUTFILES}" ENFUSE_PSTEXFILES)
add_custom_target(eps_figures_doc DEPENDS ${ENBLEND_PSTEXFILES} ${ENFUSE_PSTEXFILES})

set(INPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(OUTPUT_DIR "${CMAKE_BINARY_DIR}/doc")
configure_file(createdvi.cmake.in ${CMAKE_BINARY_DIR}/doc/createdvi.cmake @ONLY)

# build dvi/ps/pdf version of manual
add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/doc/enblend.dvi"
    COMMAND ${CMAKE_COMMAND} -Dfilename=enblend -Dcompiler=${LATEX_COMPILER} -Dflags=${LATEX_FLAGS} -P ${CMAKE_BINARY_DIR}/doc/createdvi.cmake
    DEPENDS
        ${ENBLEND_SOURCE_FILES}
        generated_files_for_doc
        eps_figures_doc
        endoc.ist
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc
)

add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/doc/enfuse.dvi"
    COMMAND ${CMAKE_COMMAND} -Dfilename=enfuse -Dcompiler=${LATEX_COMPILER} -Dflags=${LATEX_FLAGS} -P ${CMAKE_BINARY_DIR}/doc/createdvi.cmake
    DEPENDS
        ${ENFUSE_SOURCE_FILES}
        generated_files_for_doc
        eps_figures_doc
        endoc.ist
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc
)

add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/doc/enblend.ps"
    COMMAND ${DVIPS_CONVERTER} ${DVIPS_FLAGS} -o enblend.ps enblend.dvi
    DEPENDS "${CMAKE_BINARY_DIR}/doc/enblend.dvi"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc
)

add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/doc/enblend.pdf"
    COMMAND ${CMAKE_COMMAND} -Dfilename=enblend -Dcompiler=${PDFLATEX_COMPILER} -Dflags=${PDFLATEX_FLAGS} -P ${CMAKE_BINARY_DIR}/doc/createdvi.cmake
    DEPENDS
        ${ENBLEND_SOURCE_FILES}
        generated_files_for_doc
        eps_figures_doc
        endoc.ist
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc
)

add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/doc/enfuse.ps"
    COMMAND ${DVIPS_CONVERTER} ${DVIPS_FLAGS} -o enfuse.ps enfuse.dvi
    DEPENDS "${CMAKE_BINARY_DIR}/doc/enfuse.dvi"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc
)

add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/doc/enfuse.pdf"
    COMMAND ${CMAKE_COMMAND} -Dfilename=enfuse -Dcompiler=${PDFLATEX_COMPILER} -Dflags=${PDFLATEX_FLAGS} -P ${CMAKE_BINARY_DIR}/doc/createdvi.cmake
    DEPENDS
        ${ENFUSE_SOURCE_FILES}
        generated_files_for_doc
        eps_figures_doc
        endoc.ist
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc
)

# build html version of manual

# prepare picture for html version
add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/doc/seam-line-visualization.png"
    COMMAND ${IMAGEMAGICK_CONVERT_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/seam-line-visualization.tif "${CMAKE_BINARY_DIR}/doc/seam-line-visualization.png"
    DEPENDS seam-line-visualization.tif
)
set(ENBLEND_PNGFILES "${CMAKE_BINARY_DIR}/doc/seam-line-visualization.png")

macro(build_png _list _depends)
    foreach(_file ${_list})
        get_filename_component(_name ${_file} NAME_WE)
        get_filename_component(_ext ${_file} EXT)
        if(_ext STREQUAL ".svg")
            add_custom_command(
                OUTPUT "${CMAKE_BINARY_DIR}/doc/${_name}.png"
                COMMAND ${RSVGCONVERT_EXECUTABLE} ${SVGCONVERT_PNG_FLAGS} --output=${_name}.png ${CMAKE_CURRENT_SOURCE_DIR}/${_name}.svg
                DEPENDS ${_name}.svg
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc
                VERBATIM
            )
        elseif(_ext STREQUAL ".dot")
            math(EXPR _scale_png "9600 / ${PRINTRESOLUTION}")  # in makefile.am (96*1000) / print resolution / 10
            add_custom_command(
                OUTPUT "${CMAKE_BINARY_DIR}/doc/${_name}.png"
                COMMAND ${M4_EXECUTABLE} ${DOTPP_FLAGS} "--define=dot_output_type=png" ${DOTPP_PRELOAD} ${CMAKE_CURRENT_SOURCE_DIR}/${_name}.dot | ${DOT_EXECUTABLE} ${DOT_FLAGS} -Tpng ${DOT_PNG_FLAGS} -o${_name}.png
                COMMAND ${IMAGEMAGICK_CONVERT_EXECUTABLE} ${_name}.png -transparent white -resize ${_scale_png}% ${_name}.png 
                DEPENDS ${_name}.dot ${CMAKE_CURRENT_SOURCE_DIR}/uml-dot.m4 
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc
                VERBATIM
            )
            list(APPEND ${_depends} "${CMAKE_BINARY_DIR}/doc/${_name}.png")
        elseif(_ext STREQUAL ".gp")
            add_custom_command(
                OUTPUT "${CMAKE_BINARY_DIR}/doc/${_name}.png" "${CMAKE_BINARY_DIR}/doc/${_name}.cleantex.gp"
                COMMAND ${PERL_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/cleantex" ${CMAKE_CURRENT_SOURCE_DIR}/${_name}.gp > ${_name}.cleantex.gp
                COMMAND ${GNUPLOT_EXECUTABLE} ${GNUPLOT_FLAGS} -e "set output \"${_name}.png\"; set terminal pngcairo ${GNUPLOT_PNG_TERMINAL_OPTIONS}" ${_name}.cleantex.gp
                DEPENDS ${_name}.gp colors.gp
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc
                VERBATIM
            )
        else()
            # this should never happens
            message(FATAL_ERROR "Unknown extension")
        endif()
        list(APPEND ${_depends} "${CMAKE_BINARY_DIR}/doc/${_name}.png")
    endforeach()
endmacro()
build_png("${ENBLEND_IMAGE_INPUTFILES}" ENBLEND_PNGFILES)
set(ENFUSE_PNGFILES)
build_png("${ENFUSE_IMAGE_INPUTFILES}" ENFUSE_PNGFILES)
add_custom_target(png_figures_doc DEPENDS ${ENBLEND_PNGFILES} ${ENFUSE_PNGFILES})

# now build html doc
add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/doc/enblend.html"
    COMMAND ${HEVEA_EXE} ${LATEX2HTML_FLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/enblend.tex
    DEPENDS
        ${ENBLEND_SOURCE_FILES}
        generated_files_for_doc
        png_figures_doc
)
add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/doc/enfuse.html"
    COMMAND ${HEVEA_EXE} ${LATEX2HTML_FLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/enfuse.tex
    DEPENDS
        ${ENFUSE_SOURCE_FILES}
        generated_files_for_doc
        png_figures_doc
)
# hold all together
add_custom_target(ps
    DEPENDS
        "${CMAKE_BINARY_DIR}/doc/enblend.ps"
        "${CMAKE_BINARY_DIR}/doc/enfuse.ps"
)

add_custom_target(pdf
    DEPENDS
        "${CMAKE_BINARY_DIR}/doc/enblend.pdf"
        "${CMAKE_BINARY_DIR}/doc/enfuse.pdf"
)

add_custom_target(html
    DEPENDS
        "${CMAKE_BINARY_DIR}/doc/enblend.html"
        "${CMAKE_BINARY_DIR}/doc/enfuse.html"
)

# now build final rules and install targets
add_custom_target(doc ALL)

if(INSTALL_PS_DOC)
add_dependencies(doc ps)
install(FILES
        "${CMAKE_BINARY_DIR}/doc/enblend.ps"
        "${CMAKE_BINARY_DIR}/doc/enfuse.ps"
    DESTINATION ${DocumentsDirectory}
)
endif()

if(INSTALL_PDF_DOC)
add_dependencies(doc pdf)
install(FILES
        "${CMAKE_BINARY_DIR}/doc/enblend.pdf"
        "${CMAKE_BINARY_DIR}/doc/enfuse.pdf"
    DESTINATION ${DocumentsDirectory}
)
endif()

if(INSTALL_HTML_DOC)
add_dependencies(doc html)
install(FILES
        "${CMAKE_BINARY_DIR}/doc/enblend.html"
        "${CMAKE_BINARY_DIR}/doc/enfuse.html"
        ${ENBLEND_PNGFILES}
        ${ENFUSE_PNGFILES}
    DESTINATION ${DocumentsDirectory}
)
endif()
