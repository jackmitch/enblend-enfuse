info_TEXINFOS = enblend.texi \
                enfuse.texi

enblend_TEXINFOS = auxmac.texi auxmac.tex \
                   varsenblend.texi fdl.texi \
                   mask-template-characters.texi \
                   helpful-programs.texi \
                   color-profiles.texi \
                   tuning-memory-usage.texi \
                   understanding-masks.texi \
                   workflow.texi filespec.texi \
                   bug-reports.texi authors.texi \
                   photographic-workflow.fig \
                   default.css

enfuse_TEXINFOS = auxmac.texi auxmac.tex \
                  varsenfuse.texi fdl.texi \
                  mask-template-characters.texi \
                  helpful-programs.texi \
                  color-profiles.texi \
                  tuning-memory-usage.texi \
                  understanding-masks.texi \
                  workflow.texi filespec.texi \
                  bug-reports.texi authors.texi \
                  config-edge.gp config.gp \
                  entropy-cutoff.gp.in \
                  exposure-cutoff.gp.in \
                  entropy.gp.in gaussian.gp.in \
                  laplacian-of-gaussian.gp.in \
                  sharp-edge.gp.in sharp-edge.data \
                  smooth-edge.gp.in smooth-edge.data \
                  local-analysis-window.fig \
                  photographic-workflow.fig \
                  focus-stack-decision-tree.fig \
                  default.css

AM_MAKEINFOFLAGS = @AM_MAKEINFOFLAGS@ \
                   -I $(top_builddir) -I $(srcdir)
AM_MAKEINFOHTMLFLAGS = @AM_MAKEINFOHTMLFLAGS@ \
                   -I $(top_builddir) -I $(srcdir) \
                   --css-include=@srcdir@/default.css \
                   $(MAKEINFOHTMLFLAGS)
export TEXINPUTS=$(top_builddir):$(srcdir)

TEXI2DVI = texi2dvi $(TEXI2DVIFLAGS) $(EXTRATEXI2DVIFLAGS)

EXTRA_DIST = docstrings

# created by make(1), user probably wants to rebuild (often?)
MOSTLYCLEANFILES = varsenblend.texi varsenfuse.texi

# created by make(1)
CLEANFILES = entropy.{txt,eps,svg,pdf} \
             entropy-cutoff.{txt,eps,svg,pdf} \
             exposure-cutoff.{txt,eps,svg,pdf} \
             gaussian.{txt,eps,svg,pdf} \
             laplacian-of-gaussian.{txt,eps,svg,pdf} \
             local-analysis-window.{txt,eps,svg,pdf} \
             photographic-workflow.{txt,eps,svg,pdf} \
             focus-stack-decision-tree.{txt,eps,svg,pdf} \
             sharp-edge.{txt,eps,svg,pdf} \
             smooth-edge.{txt,eps,svg,pdf} \
             *.fig.bak \
             *.msg

# created by configure(1) or LaTeX
DISTCLEANFILES = texinfo.tex \
                 enblend.sc enblend.scs \
                 enfuse.sc enfuse.scs

MAINTAINERCLEANFILES = enblend.info enfuse.info


# Phony Targets

.PHONY: install-html
install-html:
	@echo '*** Target "install-html" has been deprecated.'
	@echo '*** Use target "install-xhtml" instead of "install-html".'
	false

.PHONY: html
html:
	@echo '*** Target "html" has been deprecated.'
	@echo '*** Use target "xhtml" instead of "html".'
	false

.PHONY: clean-local
clean-local:
	-rm -rf enblend.xhtml enfuse.xhtml


XMLLINT_FLAGS = --noout --valid

# The following exported variables are respected by html2xhml(1).
export TIDY SED XMLLINT

.PHONY: xhtml
xhtml: $(HTMLS:.html=.xhtml)

.PHONY: validate-xhtml
validate-xhtml: $(HTMLS:.html=.xhtml)
	for x in $(HTMLS); do \
          xhtml=`basename $$x .html`.xhtml; \
	  if test -d $$xhtml; then \
	    for y in $$xhtml/*.xhtml; do \
	      $(XMLLINT) $(XMLLINT_FLAGS) $$y; \
	    done; \
	  else \
	    $(XMLLINT) $(XMLLINT_FLAGS) $$xhtml; \
          fi; \
	done

.PHONY: install-xhtml
install-xhtml: xhtml
	@$(NORMAL_INSTALL)
	test -z "$(htmldir)" || $(MKDIR_P) "$(DESTDIR)$(htmldir)"
	for x in $(HTMLS); do \
          xhtml=`basename $$x .html`.xhtml; \
	  if test -d $$x; then \
	    $(MKDIR_P) "$(DESTDIR)$(htmldir)/$$xhtml"; \
	    $(INSTALL_DATA) $$xhtml/* "$(DESTDIR)$(htmldir)/$$xhtml"; \
	  else \
	    $(INSTALL_DATA) *.xhtml *.svg "$(DESTDIR)$(htmldir)"; \
	  fi; \
	done

.PHONY: uninstall-xhtml
uninstall-xhtml:
	@$(NORMAL_UNINSTALL)
	for x in $(HTMLS); do \
          xhtml=`basename $$x .html`.xhtml; \
	  if test -d $$x; then \
	    rm -rf "$(DESTDIR)$(htmldir)/$$xhtml"; \
	  else \
	    rm -rf "$(DESTDIR)$(htmldir)/*.xhtml" "$(DESTDIR)$(htmldir)/*.svg"; \
	  fi; \
	done


# Implicit Rules

.gp.txt:
	$(GNUPLOT) $<

.gp.pdf:
	$(GNUPLOT) $<

.gp.eps:
	$(GNUPLOT) $<

.gp.svg:
	$(GNUPLOT) $<


.fig.txt:
	$(SED) -e '1,/---BEGIN-TEXT---/d' \
               -e '/---END-TEXT---/,$$d' \
               -e 's/^# \?//' < $< > $@

.fig.eps:
	$(FIG2DEV) -L eps $< $@

.fig.pdf:
	$(FIG2DEV) -L pdf $< $@

.fig.svg:
	$(FIG2DEV) -L svg $< | \
	    $(SED) -e '/<desc>/d' -e 's/style="stroke-width:.025in/style="stroke-width:0.001px/' > $@


.html.xhtml:
	$(srcdir)/html2xhtml \
	    --document-version=$(VERSION) \
	    --tidy-flags="$(EXTRATIDYFLAGS)" \
	    --xmllint-flags="$(EXTRAXMLLINTFLAGS)" \
	    $<
	test -d "$<" && test "$<" = enblend.html && \
            for i in $(enblend_TEXINFOS); do \
              test -f $${i%%.*}.svg && cp $${i%%.*}.svg "$@"; \
            done; true
	test -d "$<" && test "$<" = enfuse.html && \
            for i in $(enfuse_TEXINFOS); do \
              test -f $${i%%.*}.svg && cp $${i%%.*}.svg "$@"; \
            done; true


# Explicit Rules

varsenblend.texi: $(srcdir)/../src/enblend.cc \
                  $(srcdir)/../src/bounds.h \
                  $(srcdir)/../src/common.h \
                  $(srcdir)/../src/global.h
	$(PERL) $(srcdir)/docstrings $^ > $@

enblend.info: $(enblend_TEXINFOS) \
              $(srcdir)/versenblend.texi $(top_builddir)/config-h.texi \
              photographic-workflow.txt

enblend.dvi: $(enblend_TEXINFOS) \
             $(srcdir)/versenblend.texi $(top_builddir)/config-h.texi \
             photographic-workflow.eps

enblend.html: $(enblend_TEXINFOS) \
              $(srcdir)/versenblend.texi $(top_builddir)/config-h.texi \
              photographic-workflow.svg

enblend.pdf: $(enblend_TEXINFOS) \
             photographic-workflow.pdf

enblend.xhtml: enblend.html


varsenfuse.texi: $(srcdir)/../src/enfuse.cc \
                 $(srcdir)/../src/bounds.h \
                 $(srcdir)/../src/common.h \
                 $(srcdir)/../src/global.h
	$(PERL) $(srcdir)/docstrings $^ > $@

enfuse.info: $(enfuse_TEXINFOS) \
             $(srcdir)/versenblend.texi $(top_builddir)/config-h.texi \
             photographic-workflow.txt \
             focus-stack-decision-tree.txt \
             entropy.txt \
             entropy-cutoff.txt \
             exposure-cutoff.txt \
             gaussian.txt \
             laplacian-of-gaussian.txt \
             local-analysis-window.txt \
             sharp-edge.txt \
             smooth-edge.txt

enfuse.dvi: $(enfuse_TEXINFOS) \
            $(srcdir)/versenfuse.texi $(top_builddir)/config-h.texi \
            photographic-workflow.eps \
            focus-stack-decision-tree.eps \
            entropy.eps \
            entropy-cutoff.eps \
            exposure-cutoff.eps \
            gaussian.eps \
            laplacian-of-gaussian.eps \
            local-analysis-window.eps \
            sharp-edge.eps \
            smooth-edge.eps

enfuse.html: $(enfuse_TEXINFOS) \
             $(srcdir)/versenfuse.texi $(top_builddir)/config-h.texi \
             photographic-workflow.svg \
             focus-stack-decision-tree.svg \
             entropy.svg \
             entropy-cutoff.svg \
             exposure-cutoff.svg \
             gaussian.svg \
             laplacian-of-gaussian.svg \
             local-analysis-window.svg \
             sharp-edge.svg \
             smooth-edge.svg

enfuse.pdf: $(enfuse_TEXINFOS) \
            $(srcdir)/versenfuse.texi $(top_builddir)/config-h.texi \
            photographic-workflow.pdf \
            focus-stack-decision-tree.pdf \
            entropy.pdf \
            entropy-cutoff.pdf \
            exposure-cutoff.pdf \
            gaussian.pdf \
            laplacian-of-gaussian.pdf \
            local-analysis-window.pdf \
            sharp-edge.pdf \
            smooth-edge.pdf

enfuse.xhtml: enfuse.html
