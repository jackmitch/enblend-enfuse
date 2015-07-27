.NOTPARALLEL:


TEXINPUTS = .:$(LATEX2HTML_LIBDIR):$(srcdir):
export TEXINPUTS

LATEX_FLAGS = -file-line-error -halt-on-error -interaction nonstopmode

## Enabled pieces of information in a LaTeX *.log file in non-verbose mode.
##    h    outputs only warnings about horizontal boxes
##    r    displays warnings about references
##    t    displays all *.tex files used
TEXLOGANALYSER_FLAGS = $(if $(subst cat,,$(TEXLOGANALYSER)),-hrt,)

MAKEINDEX_FLAGS = $(if $(subst false,,$(AM_V_P)),,-q)

## Papersize for some converters, usually:
##    a4
##    a5
##    b5
##    executive
##    legal
##    letter
PAPERSIZE = a4

## Resolution for some converters in dots-per-inch
PRINT_RESOLUTION = 600

DVIPS_FLAGS = -q -t $(PAPERSIZE) -D$(PRINT_RESOLUTION) -Z

## Some other "interesting" flags for Hevea:
##     -dv        add borders to block-level elements
##     -mathml    enable MathML output
LATEX2HTML_FLAGS = -O -I . -I $(srcdir) -fix

GNUPLOT_FLAGS = --default-settings -e 'DATA_DIR="$(srcdir)"'

## Some possible plot sizes (in pixels) with identical aspect ratios:
##     648,480
##     720,540
##     800,600
GNUPLOT_PNG_TERMINAL_OPTIONS = enhanced transparent truecolor size 720,540

## Some possible plot sizes (in inches) with identical aspect ratios:
##     4.0,3.0
##     4.5,3.375
##     5,3.75
GNUPLOT_EPSLATEX_TERMINAL_OPTIONS = color size 4.5in,3.375in

DOTPP = $(M4)
DOTPP_FLAGS = --fatal-warnings --prefix-builtins --synclines $(srcdir)/uml-dot.m4
DOT_FLAGS = \
    -Gratio=auto -Gsplines=ortho \
    -Nfontname=Helvetica -Nfontsize=11 \
    -Efontname=Helvetica -Efontsize=10
DOT_EPS_FLAGS = -Gsize=4.75,8.5
DOT_PNG_FLAGS = -Gbgcolor=transparent #-Gresolution=$(PRINT_RESOLUTION)

FIG2DEV_FLAGS =
FIG2DEV_EPS_TERMINAL_OPTIONS =
FIG2DEV_PSTEX_TERMINAL_OPTIONS =
FIG2DEV_PNG_TERMINAL_OPTIONS =

CLEANTEX = $(srcdir)/cleantex
CONFIG2TEX = $(srcdir)/config2tex
DOCSTRINGS = $(srcdir)/docstrings
UPDATED_ON = $(srcdir)/updated-on
UPDATED_ON_FLAGS = \
    --database=source-code-manager \
    --override='none:\datebyunknown' \
    --override='fs:\datebyfs' \
    --override='scm:\datebyscm'


VPATH =# must be empty to keep srcdir free of products


########################################################################


.PHONY: all-local
all-local: ps-local html-local


.PHONY: ps-local
ps-local: enblend.ps enfuse.ps


.PHONY: pdf-local
pdf-local: enblend.pdf enfuse.pdf


.PHONY: html-local
html-local: enblend.html enfuse.html


.PHONY: dvi-local
dvi-local: enblend.dvi enfuse.dvi


.PHONY: enblend
enblend: enblend.ps enblend.html


.PHONY: enfuse
enfuse: enfuse.ps enfuse.html


.PHONY: clean-local
clean-local:
	- rm -f *.eps *.png *.pstex *.gp
	- rm -f *.{g,o,p,a}nd *.{g,o,p,a}dx
	- rm -f *.idx *.ilg *.ind *.lo{a,f,t,x} *.toc
	- rm -f *.aux *.log *.dvi *.ps *.pdf
	- rm -f *.html *.haux *.hind *.hlog
	- rm -f en*.{htoc,hpnd,hond,hgnd,haux,hand,hdgx}
	- rm -f enblend.out enblend.image.{tex,out} enblend-{variables,version}.tex
	- rm -f enfuse.out enfuse.image.{tex,out} enfuse-{variables,version}.tex


.PHONY: distclean-local
distclean-local: clean-local
	- rm -f config-h.tex dynamic-preamble.tex
	- rm -f *~


.PHONY: FORCE


########################################################################


config2tex_verbose = $(config2tex_verbose_@AM_V@)
config2tex_verbose_ = $(config2tex_verbose_@AM_DEFAULT_V@)
config2tex_verbose_0 = @echo "  CONFIG2TEX $@";

convert_verbose = $(convert_verbose_@AM_V@)
convert_verbose_ = $(convert_verbose_@AM_DEFAULT_V@)
convert_verbose_0 = @echo "  CONVERT  $@";

docstrings_verbose = $(docstrings_verbose_@AM_V@)
docstrings_verbose_ = $(docstrings_verbose_@AM_DEFAULT_V@)
docstrings_verbose_0 = @echo "  DOCSTR   $@";

dvips_verbose = $(dvips_verbose_@AM_V@)
dvips_verbose_ = $(dvips_verbose_@AM_DEFAULT_V@)
dvips_verbose_0 = @echo "  DVIPS    $@";

dot_verbose = $(dot_verbose_@AM_V@)
dot_verbose_ = $(dot_verbose_@AM_DEFAULT_V@)
dot_verbose_0 = @echo "  DOT      $@";

fig2dev_verbose = $(fig2dev_verbose_@AM_V@)
fig2dev_verbose_ = $(fig2dev_verbose_@AM_DEFAULT_V@)
fig2dev_verbose_0 = @echo "  FIG2DEV  $@";

gnuplot_verbose = $(gnuplot_verbose_@AM_V@)
gnuplot_verbose_ = $(gnuplot_verbose_@AM_DEFAULT_V@)
gnuplot_verbose_0 = @echo "  GNUPLOT  $@";

install_verbose = $(install_verbose_@AM_V@)
install_verbose_ = $(install_verbose_@AM_DEFAULT_V@)
install_verbose_0 = @echo "  INSTALL  $@";

latex_verbose = $(call latex_verbose_@AM_V@,$(1))
latex_verbose_ = $(call latex_verbose_@AM_DEFAULT_V@,$(1))
latex_verbose_0 = @echo "  LATEX    [pass $(1)] $@";

makeindex_verbose = $(call makeindex_verbose_@AM_V@,$(1))
makeindex_verbose_ = $(call makeindex_verbose_@AM_DEFAULT_V@,$(1))
makeindex_verbose_0 = @echo "  MKIDX    [$(1)] $@";


########################################################################


## Alternative PDF generation: LaTeX -> DVI -> PostScript -> PDF
PS2PDF_FLAGS = -r$(PRINT_RESOLUTION)
%.ps.pdf: %.ps
	$(AM_V_GEN)$(PS2PDF) $(PS2PDF_FLAGS) $< $@


## Alternative PDF generation: LaTeX -> DVI -> PDF
DVIPDF_FLAGS = -p $(PAPERSIZE) -r $(PRINT_RESOLUTION) -v
%.dvi.pdf: %.dvi
	$(AM_V_GEN)$(DVIPDF) $(DVIPDF_FLAGS) -o $@ $<


%.ps: %.dvi
	$(dvips_verbose)$(DVIPS) $(DVIPS_FLAGS) -o $@ $<


SHOW_LATEX_LOG = \
    rc=$$?; \
    if $(AM_V_P) || test $$rc -ne 0; \
    then \
      cat $(1); \
    else \
      $(TEXLOGANALYSER) $(TEXLOGANALYSER_FLAGS) $(1); \
    fi; \
    test $$rc -eq 0

%.dvi %.pdf: $(srcdir)/%.tex $(srcdir)/endoc.ist
	@- rm -f $(foreach ext,aux toc lo[^g] ?nd ?dx,$(basename $@).$(ext))
	$(call latex_verbose,1)$(LATEX_TRANSLATOR) $(LATEX_FLAGS) $(EXTRA_LATEX_FLAGS) $< >/dev/null; \
            $(call SHOW_LATEX_LOG,$*.log)
	$(call makeindex_verbose,general terms)$(MAKEINDEX) $(MAKEINDEX_FLAGS) -s $(srcdir)/endoc.ist -o $*.gnd $*.gdx
	$(call makeindex_verbose,command-line options)$(MAKEINDEX) $(MAKEINDEX_FLAGS) -o $*.ond $*.odx
	$(call makeindex_verbose,applications and programs)$(MAKEINDEX) $(MAKEINDEX_FLAGS) -o $*.pnd $*.pdx
	$(call makeindex_verbose,syntactic comments)$(MAKEINDEX) $(MAKEINDEX_FLAGS) -o $*.and $*.adx
	$(call latex_verbose,2)$(LATEX_TRANSLATOR) $(LATEX_FLAGS) $(EXTRA_LATEX_FLAGS) $< >/dev/null; \
            $(call SHOW_LATEX_LOG,$*.log)
	@ max_iter=5; \
          iter=2; \
          while [ $$iter -le $$max_iter ] && \
                $(GREP) '^LaTeX Warning: There were undefined references' $*.log >/dev/null && \
                $(GREP) '^LaTeX Warning: Label(s) may have changed' $*.log >/dev/null; \
          do \
            iter=$$((iter + 1)); \
            if $(AM_V_P); then :; else printf '  LATEX    [pass %d] %s\n' $$iter $@; fi; \
            $(LATEX_TRANSLATOR) $(LATEX_FLAGS) $(EXTRA_LATEX_FLAGS) $< >/dev/null; \
                $(call SHOW_LATEX_LOG,$*.log); \
          done; \
          if [ $$iter -le $$max_iter ]; \
          then \
            printf 'Fixpoint reached in %i step(s)\n' $$iter; \
          else \
            printf 'Fixpoint NOT reached after %i steps; still undefined or loose labels\n' $$iter; \
            $(GREP) '^LaTeX Warning: Reference' $*.log; \
            false; \
          fi


SHOW_LATEX2HTML_LOG = \
    rc=$$?; \
    if $(AM_V_P); \
    then \
      cat $(1); \
    else \
      $(EGREP) '(Warning:|^Fixpoint)' $(1); \
    fi; \
    test $$rc -eq 0

%.html: $(srcdir)/%.tex
	$(AM_V_GEN)$(LATEX2HTML) $(LATEX2HTML_FLAGS) $< >$*.hlog 2>&1; \
            $(call SHOW_LATEX2HTML_LOG,$*.hlog)


%.cleantex.gp: $(srcdir)/%.gp
	$(AM_V_at)$(CLEANTEX) $^ > $@


%.png: %.cleantex.gp
	$(gnuplot_verbose)$(GNUPLOT) $(GNUPLOT_FLAGS) \
            -e 'set output "$@"; set terminal pngcairo $(GNUPLOT_PNG_TERMINAL_OPTIONS)' $<


%.pstex: $(srcdir)/%.gp
	$(gnuplot_verbose)$(GNUPLOT) $(GNUPLOT_FLAGS) \
            -e 'set output "$*.tex"; set terminal epslatex $(GNUPLOT_EPSLATEX_TERMINAL_OPTIONS)' $< && \
        mv $*.tex $@


%.eps: $(srcdir)/%.dot
	$(dot_verbose)$(DOTPP) $(DOTPP_FLAGS) $< | $(DOT) $(DOT_FLAGS) -Teps $(DOT_EPS_FLAGS) -o$@

%.png: $(srcdir)/%.dot
	$(dot_verbose)$(DOTPP) $(DOTPP_FLAGS) $< | $(DOT) $(DOT_FLAGS) -Tpng $(DOT_PNG_FLAGS) -o$@


%.eps: $(srcdir)/%.fig
	$(fig2dev_verbose)$(FIG2DEV) $(FIG2DEV_FLAGS) -L eps $(FIG2DEV_EPS_TERMINAL_OPTIONS) $< $@


%.pstex: $(srcdir)/%.fig
	$(fig2dev_verbose)$(FIG2DEV) $(FIG2DEV_FLAGS) -L pstex $(FIG2DEV_EPS_TERMINAL_OPTIONS) $< $*.eps
	$(fig2dev_verbose)$(FIG2DEV) $(FIG2DEV_FLAGS) -L pstex_t -p $*.eps $(FIG2DEV_PSTEX_TERMINAL_OPTIONS) $< $@


# We save at 4x scale to downscale with convert's fine aliasing capabilities.
# Order of convert options is important: First throw away the white
# part as transparent background and afterwards resize with
# intermediate gray-values.
%.png: $(srcdir)/%.fig
	$(fig2dev_verbose)$(FIG2DEV) $(FIG2DEV_FLAGS) -L png $(FIG2DEV_PNG_TERMINAL_OPTIONS) -m 4 $< | \
            $(CONVERT) png:- -transparent white -resize 30% $@


%.eps: $(srcdir)/%.tif
	$(convert_verbose)$(CONVERT) $< eps3:$@


%.png: $(srcdir)/%.tif
	$(convert_verbose)$(CONVERT) $< $@


%-version.tex: FORCE
	$(docstrings_verbose)$(DOCSTRINGS) key-value \
            UPDATED '$(shell $(UPDATED_ON) $(UPDATED_ON_FLAGS) $(FILES))' \
            VERSION '$(VERSION)' \
            DATAROOTDIR '$(datarootdir)' > $@


########################################################################


DOCSTRING_SOURCE_FILES = \
    $(addprefix $(srcdir)/../src/, \
        bounds.h \
        common.h \
        filespec.h filespec.cc \
        global.h \
        mga.h \
        opencl.cc \
        layer_selection/selector.cc \
    )


config-h.tex: ../config.h
	$(config2tex_verbose)$(CONFIG2TEX) $< $@


dynamic-preamble.tex: FORCE
	$(AM_V_GEN)
	$(AM_V_at)printf '%%%%%% %s\n' '$(basename $@)' > $@
	$(AM_V_at)printf '\PassOptionsToPackage{%spaper}{report}\n' '$(PAPERSIZE)' >> $@
	$(AM_V_at)printf '\PassOptionsToPackage{%spaper}{refrep}\n' '$(PAPERSIZE)' >> $@
	$(AM_V_at)printf '\\newif\ifreferencemanual\n' >> $@
	$(AM_V_at)printf '\\newif\ifhyperref\n' >> $@
	$(AM_V_at)printf '%%%% User definitions follow...\n' >> $@
	$(AM_V_at)printf '%s\n' '$(DYNAMIC_TEX_PREAMBLE)' >> $@


external-mask-workflow.eps focus-stack-decision-tree.eps photographic-workflow.eps \
external-mask-workflow.png focus-stack-decision-tree.png photographic-workflow.png: $(srcdir)/uml-dot.m4

entropy-cutoff.gp entropy.gp exposure-cutoff.gp exposure-weights.gp \
fullsine.gp gaussian.gp halfsine.gp laplacian-of-gaussian.gp \
log-transform.gp lorentzian.gp power.gp: colors.gp


########################################################################


ENBLEND_SOURCE_FILES = \
    $(wildcard $(srcdir)/common-*.tex) \
    $(wildcard $(srcdir)/enblend-*.tex) \
    $(addprefix $(srcdir)/, \
        enblend.tex \
        static-preamble.tex  lead-in.tex  lead-out.tex \
        index-def.tex  index-use.tex \
        macros.tex  exemplar.hva \
        htmlstyle.tex  url-def.tex \
    )


ENBLEND_HYBRID_MATERIAL = \
    log-transform.gp


ENBLEND_EXTRA_MATERIAL = \
    external-mask-workflow.dot \
    photographic-workflow.dot \
    \
    seam-line-visualization.tif


ENBLEND_GENERATED_FILES = \
    config-h.tex \
    dynamic-preamble.tex \
    enblend-variables.tex \
    enblend-version.tex


ENBLEND_DVI_SEPARATE_PREREQUISITES = \
    $(addsuffix .pstex,$(basename $(ENBLEND_HYBRID_MATERIAL))) \
    $(addsuffix .eps,$(basename $(ENBLEND_EXTRA_MATERIAL)))


enblend.pdf: LATEX_TRANSLATOR = $(PDFLATEX)

enblend.dvi: LATEX_TRANSLATOR = $(LATEX)

enblend.pdf enblend.dvi: \
    $(ENBLEND_SOURCE_FILES) \
    $(ENBLEND_GENERATED_FILES) \
    $(ENBLEND_DVI_SEPARATE_PREREQUISITES)


ENBLEND_HTML_SEPARATE_PREREQUISITES = \
    $(addsuffix .png,$(basename $(ENBLEND_HYBRID_MATERIAL))) \
    $(addsuffix .png,$(basename $(ENBLEND_EXTRA_MATERIAL)))


enblend.html: \
    $(ENBLEND_SOURCE_FILES) \
    $(ENBLEND_GENERATED_FILES) \
    $(ENBLEND_HTML_SEPARATE_PREREQUISITES)


enblend-variables.tex: $(DOCSTRING_SOURCE_FILES) $(srcdir)/../src/enblend.cc
	$(docstrings_verbose)$(DOCSTRINGS) scan-files $^ > $@


enblend-version.tex: FILES = $(ENBLEND_SOURCE_FILES)


########################################################################


ENFUSE_SOURCE_FILES = \
    $(wildcard $(srcdir)/common-*.tex) \
    $(wildcard $(srcdir)/enfuse-*.tex) \
    $(addprefix $(srcdir)/, \
        enfuse.tex \
        static-preamble.tex  lead-in.tex  lead-out.tex \
        index-def.tex  index-use.tex \
        macros.tex  exemplar.hva \
        htmlstyle.tex  url-def.tex \
    )

ENFUSE_HYBRID_MATERIAL = \
    local-analysis-window.fig \
    \
    entropy-cutoff.gp \
    entropy.gp \
    exposure-cutoff.gp \
    exposure-weights.gp \
    fullsine.gp \
    gaussian.gp \
    halfsine.gp \
    laplacian-of-gaussian.gp \
    log-transform.gp \
    lorentzian.gp \
    power.gp \
    sharp-edge.gp \
    smooth-edge.gp


ENFUSE_EXTRA_MATERIAL = \
    external-mask-workflow.dot \
    focus-stack-decision-tree.dot \
    photographic-workflow.dot


ENFUSE_GENERATED_FILES = \
    config-h.tex \
    dynamic-preamble.tex \
    enfuse-variables.tex \
    enfuse-version.tex


ENFUSE_DVI_SEPARATE_PREREQUISITES = \
    $(addsuffix .pstex,$(basename $(ENFUSE_HYBRID_MATERIAL))) \
    $(addsuffix .eps,$(basename $(ENFUSE_EXTRA_MATERIAL)))


enfuse.pdf: LATEX_TRANSLATOR = $(PDFLATEX)

enfuse.dvi: LATEX_TRANSLATOR = $(LATEX)

enfuse.pdf enfuse.dvi: \
    $(ENFUSE_SOURCE_FILES) \
    $(ENFUSE_GENERATED_FILES) \
    $(ENFUSE_DVI_SEPARATE_PREREQUISITES)


ENFUSE_HTML_SEPARATE_PREREQUISITES = \
    $(addsuffix .png,$(basename $(ENFUSE_HYBRID_MATERIAL))) \
    $(addsuffix .png,$(basename $(ENFUSE_EXTRA_MATERIAL)))


enfuse.html: \
    $(ENFUSE_SOURCE_FILES) \
    $(ENFUSE_GENERATED_FILES) \
    $(ENFUSE_HTML_SEPARATE_PREREQUISITES)


enfuse-variables.tex: $(DOCSTRING_SOURCE_FILES) $(srcdir)/../src/enfuse.cc
	$(docstrings_verbose)$(DOCSTRINGS) scan-files $^ > $@


enfuse-version.tex: FILES = $(ENFUSE_SOURCE_FILES)


########################################################################


.PHONY: install-ps-local
install-ps-local: ps-local
	@ $(MKDIR_P) $(psdir)
	$(install_verbose)$(INSTALL_DATA) enblend.ps enfuse.ps $(psdir)


.PHONY: install-html-local
install-html-local: html-local
	@ $(MKDIR_P) $(htmldir)
	$(install_verbose)$(INSTALL_DATA) \
            $(sort \
                enblend.html $(ENBLEND_HTML_SEPARATE_PREREQUISITES) \
                enfuse.html $(ENFUSE_HTML_SEPARATE_PREREQUISITES) \
            ) \
            $(htmldir)


.PHONY: install-pdf-local
install-pdf-local: pdf-local
	@ $(MKDIR_P) $(pdfdir)
	$(install_verbose)$(INSTALL_DATA) enblend.pdf enfuse.pdf $(pdfdir)


.PHONY: install-dvi-local
install-dvi-local: dvi-local
	@ $(MKDIR_P) $(dvidir)
	$(install_verbose)$(INSTALL_DATA) \
            $(sort \
                enblend.dvi $(patsubst %.pstex,%.eps,$(ENBLEND_DVI_SEPARATE_PREREQUISITES)) \
                enfuse.dvi $(patsubst %.pstex,%.eps,$(ENFUSE_DVI_SEPARATE_PREREQUISITES)) \
            ) \
            $(dvidir)


EXTRA_DIST = \
    $(wildcard $(srcdir)/*.data) \
    $(wildcard $(srcdir)/*.dot) \
    $(wildcard $(srcdir)/*.fig) \
    $(wildcard $(srcdir)/*.gp) \
    $(wildcard $(srcdir)/*.pm) \
    $(wildcard $(srcdir)/*.tex) \
    $(wildcard $(srcdir)/*.tif) \
    $(addprefix $(srcdir)/, \
        uml-dot.m4 \
        endoc.ist exemplar.hva \
        cleantex config2tex docstrings updated-on \
    )
