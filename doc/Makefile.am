TEXINPUTS = .:$(LATEX2HTML_LIBDIR):$(srcdir):
export TEXINPUTS

LATEX_FLAGS = -file-line-error -halt-on-error -interaction nonstopmode

MAKEINDEX_FLAGS =

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
    --override='scm:\datebyscm' \
    --format='%B %-m, %Y%q'


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


## Alternative PDF generation: LaTeX -> DVI -> PostScript -> PDF
PS2PDF_FLAGS = -r$(PRINT_RESOLUTION)
%.ps.pdf: %.ps
	$(PS2PDF) $(PS2PDF_FLAGS) $< $@


## Alternative PDF generation: LaTeX -> DVI -> PDF
DVIPDF_FLAGS = -p $(PAPERSIZE) -r $(PRINT_RESOLUTION) -v
%.dvi.pdf: %.dvi
	$(DVIPDF) $(DVIPDF_FLAGS) -o $@ $<


%.ps: %.dvi
	$(DVIPS) $(DVIPS_FLAGS) -o $@ $<


%.dvi %.pdf: $(srcdir)/%.tex $(srcdir)/endoc.ist
	- rm -f $(foreach ext,aux toc lo[^g] ?nd ?dx,$(basename $@).$(ext))
	$(LATEX_TRANSLATOR) $(LATEX_FLAGS) $(EXTRA_LATEX_FLAGS) $<
	$(MAKEINDEX) $(MAKEINDEX_FLAGS) -s $(srcdir)/endoc.ist -o $*.gnd $*.gdx
	$(MAKEINDEX) $(MAKEINDEX_FLAGS) -o $*.ond $*.odx
	$(MAKEINDEX) $(MAKEINDEX_FLAGS) -o $*.pnd $*.pdx
	$(MAKEINDEX) $(MAKEINDEX_FLAGS) -o $*.and $*.adx
	$(LATEX_TRANSLATOR) $(LATEX_FLAGS) $(EXTRA_LATEX_FLAGS) $<
	@ max_iter=5; \
	  iter=2; \
	  while [ $$iter -le $$max_iter ] && \
                grep --quiet '^LaTeX Warning: There were undefined references' $*.log && \
                grep --quiet '^LaTeX Warning: Label(s) may have changed' $*.log; \
          do \
	    $(LATEX_TRANSLATOR) $(LATEX_FLAGS) $(EXTRA_LATEX_FLAGS) $<; \
	    iter=$$((iter + 1)); \
	  done; \
	  if [ $$iter -le $$max_iter ]; \
          then \
	    printf 'Fixpoint reached in %i step(s)\n' $$iter; \
	  else \
	    printf 'Fixpoint NOT reached after %i steps; still undefined or loose labels\n' $$iter; \
	    grep '^LaTeX Warning: Reference' $*.log; \
	    false; \
	  fi


%.html: $(srcdir)/%.tex
	$(LATEX2HTML) $(LATEX2HTML_FLAGS) $< 2>&1 | tee $*.hlog


%.cleantex.gp: $(srcdir)/%.gp
	$(CLEANTEX) $^ > $@


%.png: %.cleantex.gp
	$(GNUPLOT) $(GNUPLOT_FLAGS) -e 'set output "$@"; set terminal pngcairo $(GNUPLOT_PNG_TERMINAL_OPTIONS)' $<


%.pstex: $(srcdir)/%.gp
	$(GNUPLOT) $(GNUPLOT_FLAGS) -e 'set output "$*.tex"; set terminal epslatex $(GNUPLOT_EPSLATEX_TERMINAL_OPTIONS)' $< && mv $*.tex $@


%.eps: $(srcdir)/%.fig
	$(FIG2DEV) $(FIG2DEV_FLAGS) -L eps $(FIG2DEV_EPS_TERMINAL_OPTIONS) $< $@


%.pstex: $(srcdir)/%.fig
	$(FIG2DEV) $(FIG2DEV_FLAGS) -L pstex $(FIG2DEV_EPS_TERMINAL_OPTIONS) $< $*.eps
	$(FIG2DEV) $(FIG2DEV_FLAGS) -L pstex_t -p $*.eps $(FIG2DEV_PSTEX_TERMINAL_OPTIONS) $< $@


# We save at 4x scale to downscale with convert's fine aliasing capabilities.
# Order of convert options is important: First throw away the white
# part as transparent background and afterwards resize with
# intermediate gray-values.
%.png: $(srcdir)/%.fig
	$(FIG2DEV) $(FIG2DEV_FLAGS) -L png $(FIG2DEV_PNG_TERMINAL_OPTIONS) -m 4 $< | \
	    $(CONVERT) png:- -transparent white -resize 30% $@


%.eps: $(srcdir)/%.tif
	$(CONVERT) $< eps3:$@


%.png: $(srcdir)/%.tif
	$(CONVERT) $< $@


%-version.tex: FORCE
	$(DOCSTRINGS) key-value \
	    UPDATED '$(shell $(UPDATED_ON) $(UPDATED_ON_FLAGS) $(FILES))' \
	    VERSION '$(VERSION)' \
	    DATAROOTDIR '$(datarootdir)' > $@


########################################################################


DOCSTRING_SOURCE_FILES = \
    $(addprefix $(srcdir)/../src/, \
        bounds.h \
        common.h \
        filespec.cc \
        global.h \
        mga.h \
        opencl.cc \
    )


config-h.tex: ../config.h
	$(CONFIG2TEX) $< $@


dynamic-preamble.tex: FORCE
	printf '%%%%%% %s\n' '$(basename $@)' > $@
	printf '\PassOptionsToPackage{%spaper}{report}\n' '$(PAPERSIZE)' >> $@
	printf '\PassOptionsToPackage{%spaper}{refrep}\n' '$(PAPERSIZE)' >> $@
	printf '\\newif\ifreferencemanual\n' >> $@
	printf '\\newif\ifhyperref\n' >> $@
	printf '%%%% User definitions follow...\n' >> $@
	printf '%s\n' '$(DYNAMIC_TEX_PREAMBLE)' >> $@


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
    external-mask-workflow.fig \
    photographic-workflow.fig \
    \
    log-transform.gp


ENBLEND_EXTRA_MATERIAL = \
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
	$(DOCSTRINGS) scan-files $^ > $@


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
    external-mask-workflow.fig \
    focus-stack-decision-tree.fig \
    local-analysis-window.fig \
    photographic-workflow.fig \
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


ENFUSE_EXTRA_MATERIAL =


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
	$(DOCSTRINGS) scan-files $^ > $@


enfuse-version.tex: FILES = $(ENFUSE_SOURCE_FILES)


########################################################################


.PHONY: install-ps-local
install-ps-local: ps-local
	$(MKDIR_P) $(psdir)
	$(INSTALL_DATA) enblend.ps enfuse.ps $(psdir)


.PHONY: install-html-local
install-html-local: html-local
	$(MKDIR_P) $(htmldir)
	$(INSTALL_DATA) \
	    $(sort \
	        enblend.html $(ENBLEND_HTML_SEPARATE_PREREQUISITES) \
	        enfuse.html $(ENFUSE_HTML_SEPARATE_PREREQUISITES) \
            ) \
	    $(htmldir)


.PHONY: install-pdf-local
install-pdf-local: pdf-local
	$(MKDIR_P) $(pdfdir)
	$(INSTALL_DATA) enblend.pdf enfuse.pdf $(pdfdir)


.PHONY: install-dvi-local
install-dvi-local: dvi-local
	$(MKDIR_P) $(dvidir)
	$(INSTALL_DATA) \
	    $(sort \
	    	enblend.dvi $(patsubst %.pstex,%.eps,$(ENBLEND_DVI_SEPARATE_PREREQUISITES)) \
	    	enfuse.dvi $(patsubst %.pstex,%.eps,$(ENFUSE_DVI_SEPARATE_PREREQUISITES)) \
            ) \
	    $(dvidir)


EXTRA_DIST = \
    $(wildcard $(srcdir)/*.data) \
    $(wildcard $(srcdir)/*.fig) \
    $(wildcard $(srcdir)/*.gp) \
    $(wildcard $(srcdir)/*.pm) \
    $(wildcard $(srcdir)/*.tex) \
    $(wildcard $(srcdir)/*.tif) \
    $(addprefix $(srcdir)/, \
        endoc.ist exemplar.hva \
        cleantex config2tex docstrings updated-on \
    )
