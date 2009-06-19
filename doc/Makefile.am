info_TEXINFOS = enblend.texi \
                enfuse.texi

enblend_TEXINFOS = auxmac.texi auxmac.tex \
                   fdl.texi \
                   helpful-programs.texi \
                   tuning-memory-usage.texi \
                   understanding-masks.texi \
                   workflow.texi \
                   authors.texi \
                   photographic-workflow.fig \
                   default.css

enfuse_TEXINFOS = auxmac.texi auxmac.tex \
                  fdl.texi \
                  helpful-programs.texi \
                  tuning-memory-usage.texi \
                  understanding-masks.texi \
                  workflow.texi \
                  authors.texi \
                  config-edge.gp config.gp \
                  entropy-cutoff.gp.in \
                  entropy.gp.in gaussian.gp.in \
                  laplacian-of-gaussian.gp.in \
                  sharp-edge.gp.in sharp-edge.data \
                  smooth-edge.gp.in smooth-edge.data \
                  local-analysis-window.fig \
                  photographic-workflow.fig \
                  default.css

AM_MAKEINFOFLAGS = @AM_MAKEINFOFLAGS@ \
                   --css-include=@srcdir@/default.css

# created by make(1), user pobably wants to rebuild (often?)
MOSTLYCLEANFILES =

# created by make(1)
CLEANFILES = entropy.{txt,png,eps,pdf} \
             entropy-cutoff.{txt,png,eps,pdf} \
             gaussian.{txt,png,eps,pdf} \
             laplacian-of-gaussian.{txt,png,eps,pdf} \
             local-analysis-window.{txt,png,eps,pdf} \
             photographic-workflow.{txt,png,eps,pdf} \
             sharp-edge.{txt,png,eps,pdf} \
             smooth-edge.{txt,png,eps,pdf} \
             *.fig.bak \
             *.xhtml *.msg

# created by configure(1)
DISTCLEANFILES = texinfo.tex

MAINTAINERCLEANFILES = enblend.info enfuse.info


# Phony Targets

# The public id is in "tidy.cfg".
DTD_SYSID="http://www.w3.org/TR/MathML2/dtd/xhtml-math11-f.dtd"

.PHONY: xhtml
xhtml: html
	for x in *.html; do \
	  if test -d $$x; then \
	    for y in $$x/*.html; do \
	      tidy -config @srcdir@/tidy.cfg -file $${y/%.html/.msg} $$y | \
	        $(SED) --expression='1,9s#""#$(DTD_SYSID)#' > $${y/%.html/.xhtml}; \
            done; \
	  else \
	    tidy -config @srcdir@/tidy.cfg -file $${x/%.html/.msg} $$x | \
	      $(SED) --expression='1,9s#""#$(DTD_SYSID)#' > $${x/%.html/.xhtml}; \
          fi; \
	done


# Implicit Rules

.gp.txt:
	GDFONTPATH=$${GDFONTPATH:-@GDFONTPATH@} $(GNUPLOT) $<

.gp.png:
	GDFONTPATH=$${GDFONTPATH:-@GDFONTPATH@} $(GNUPLOT) $<

.gp.eps:
	GDFONTPATH=$${GDFONTPATH:-@GDFONTPATH@} $(GNUPLOT) $<

.gp.pdf:
	GDFONTPATH=$${GDFONTPATH:-@GDFONTPATH@} $(GNUPLOT) $<


.fig.txt:
	$(SED) --expression='1,/---BEGIN-TEXT---/d' \
               --expression='/---END-TEXT---/,$$d' \
               --expression='s/^# \?//' < $< > $@

.fig.png:
	$(FIG2DEV) -L png $< $@

.fig.eps:
	$(FIG2DEV) -L eps $< $@

.fig.pdf:
	$(FIG2DEV) -L pdf $< $@


# Explicit Rules

enblend.info: photographic-workflow.txt

enblend.dvi: photographic-workflow.eps

enblend.html \
enblend.pdf: photographic-workflow.png


enfuse.info: photographic-workflow.txt \
             entropy.txt \
             entropy-cutoff.txt \
             gaussian.txt \
             laplacian-of-gaussian.txt \
             local-analysis-window.txt \
             sharp-edge.txt \
             smooth-edge.txt

enfuse.dvi: photographic-workflow.eps \
            entropy.eps \
            entropy-cutoff.eps \
            gaussian.eps \
            laplacian-of-gaussian.eps \
            local-analysis-window.eps \
            sharp-edge.eps \
            smooth-edge.eps

enfuse.html \
enfuse.pdf: photographic-workflow.png \
            entropy.png \
            entropy-cutoff.png \
            gaussian.png \
            laplacian-of-gaussian.png \
            local-analysis-window.png \
            sharp-edge.png \
            smooth-edge.png
