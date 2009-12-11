info_TEXINFOS = enblend.texi \
                enfuse.texi

enblend_TEXINFOS = auxmac.texi auxmac.tex \
                   varsenblend.texi fdl.texi \
                   mask-template-characters.texi \
                   helpful-programs.texi \
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
                  tuning-memory-usage.texi \
                  understanding-masks.texi \
                  workflow.texi filespec.texi \
                  bug-reports.texi authors.texi \
                  config-edge.gp config.gp \
                  entropy-cutoff.gp.in \
                  entropy.gp.in gaussian.gp.in \
                  laplacian-of-gaussian.gp.in \
                  sharp-edge.gp.in sharp-edge.data \
                  smooth-edge.gp.in smooth-edge.data \
                  local-analysis-window.fig \
                  photographic-workflow.fig \
                  focus-stack-decision-tree.fig \
                  default.css

AM_MAKEINFOFLAGS = @AM_MAKEINFOFLAGS@ \
                   -I $(top_builddir) -I $(srcdir) \
                   --css-include=@srcdir@/default.css
export TEXINPUTS=$(top_builddir):$(srcdir)

EXTRA_DIST = docstrings

# created by make(1), user probably wants to rebuild (often?)
MOSTLYCLEANFILES = varsenblend.texi varsenfuse.texi

# created by make(1)
CLEANFILES = entropy.{txt,png,eps,pdf} \
             entropy-cutoff.{txt,png,eps,pdf} \
             gaussian.{txt,png,eps,pdf} \
             laplacian-of-gaussian.{txt,png,eps,pdf} \
             local-analysis-window.{txt,png,eps,pdf} \
             photographic-workflow.{txt,png,eps,pdf} \
             focus-stack-decision-tree.{txt,png,eps,pdf} \
             sharp-edge.{txt,png,eps,pdf} \
             smooth-edge.{txt,png,eps,pdf} \
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

.PHONY: html-local
html-local:
	test -d enblend.html && \
            for i in $(enblend_TEXINFOS); do \
              test -f $${i%%.*}.png && cp $${i%%.*}.png enblend.html; \
            done; true
	test -d enfuse.html && \
            for i in $(enfuse_TEXINFOS); do \
              test -f $${i%%.*}.png && cp $${i%%.*}.png enfuse.html; \
            done; true

.PHONY: clean-local
clean-local:
	-rm -rf enblend.xhtml enfuse.xhtml


TIDY_FLAGS=-config $(srcdir)/tidy.cfg -file /dev/null
SED_FLAGS_XHTML=-f $(srcdir)/tidy-post.sed

.PHONY: xhtml
xhtml: $(HTMLS) html-local
	for x in $(HTMLS); do \
          xhtml=$${x/%.html/.xhtml}; \
	  if test -d $$x; then \
            test -d $$xhtml || mkdir $$xhtml; \
            cp *.png $$xhtml; \
	    for y in $$x/*.html; do \
              if test "$$x/" = "$${y%Table_*}" -o "$$x/" = "$${y%Figure_*}"; then \
                  $(SED) -e 's|>| />|' < $$y > $${y//.html/.xhtml}; \
              else \
		  $(TIDY) $(TIDY_FLAGS) $$y | \
		    $(SED) $(SED_FLAGS_XHTML) \
			   -e 's|^<head>|<head>\n<meta content="$(VERSION)" name="enblend-version" />|' \
			   -e 's|<a\([^>]*\)name="[^"]*"\([^>]*\)|<a\1\2|g' \
			   -e 's|\("[^/"]*\)\.html\([#"]\)|\1.xhtml\2|' \
			   > $${y//.html/.xhtml}; \
              fi \
            done; \
	  else \
	    $(TIDY) $(TIDY_FLAGS) $$x | \
	      $(SED) $(SED_FLAGS_XHTML) \
	    	 -e 's|^<head>|<head>\n<meta content="$(VERSION)" name="enblend-version" />|' \
	    	 -e 's|<a\([^>]*\)name=\([^>]*\)|<a\1id=\2|g' > $$xhtml; \
          fi; \
	done

.PHONY: install-xhtml
install-xhtml: xhtml
	@$(NORMAL_INSTALL)
	test -z "$(htmldir)" || $(MKDIR_P) "$(DESTDIR)$(htmldir)"
	for x in $(HTMLS); do \
          xhtml=$${x/%.html/.xhtml}; \
	  if test -d $$x; then \
	    test -d "$(DESTDIR)$(htmldir)/$$xhtml" || mkdir "$(DESTDIR)$(htmldir)/$$xhtml"; \
	    $(INSTALL_DATA) $$xhtml/* "$(DESTDIR)$(htmldir)/$$xhtml"; \
	  else \
	    $(INSTALL_DATA) *.xhtml *.png "$(DESTDIR)$(htmldir)"; \
	  fi; \
	done

.PHONY: uninstall-xhtml
uninstall-xhtml:
	@$(NORMAL_UNINSTALL)
	for x in $(HTMLS); do \
          xhtml=$${x/%.html/.xhtml}; \
	  if test -d $$x; then \
	    test -d "$(DESTDIR)$(htmldir)/$$xhtml" || mkdir "$(DESTDIR)$(htmldir)/$$xhtml"; \
	    echo rm -rf "$(DESTDIR)$(htmldir)/$$xhtml"; \
	  else \
	    echo rm -rf "$(DESTDIR)$(htmldir)/*.xhtml" "$(DESTDIR)$(htmldir)/*.png"; \
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
	$(SED) -e '1,/---BEGIN-TEXT---/d' \
               -e '/---END-TEXT---/,$$d' \
               -e 's/^# \?//' < $< > $@

.fig.png:
	$(FIG2DEV) -L png $< $@

.fig.eps:
	$(FIG2DEV) -L eps $< $@

.fig.pdf:
	$(FIG2DEV) -L pdf $< $@


# Explicit Rules

varsenblend.texi: $(srcdir)/../src/enblend.cc \
                  $(srcdir)/../src/bounds.h \
                  $(srcdir)/../src/common.h \
                  $(srcdir)/../src/global.h
	$(PERL) $(srcdir)/docstrings $^ > $@

enblend.info: photographic-workflow.txt

enblend.dvi: enblend.texi \
             $(srcdir)/versenblend.texi $(top_builddir)/config-h.texi \
             $(enblend_TEXINFOS) \
             photographic-workflow.eps

enblend.html \
enblend.pdf: photographic-workflow.png


varsenfuse.texi: $(srcdir)/../src/enfuse.cc \
                 $(srcdir)/../src/bounds.h \
                 $(srcdir)/../src/common.h \
                 $(srcdir)/../src/global.h
	$(PERL) $(srcdir)/docstrings $^ > $@

enfuse.info: photographic-workflow.txt \
             focus-stack-decision-tree.txt \
             entropy.txt \
             entropy-cutoff.txt \
             gaussian.txt \
             laplacian-of-gaussian.txt \
             local-analysis-window.txt \
             sharp-edge.txt \
             smooth-edge.txt

enfuse.dvi: enfuse.texi \
            photographic-workflow.eps \
            focus-stack-decision-tree.eps \
            entropy.eps \
            entropy-cutoff.eps \
            gaussian.eps \
            laplacian-of-gaussian.eps \
            local-analysis-window.eps \
            sharp-edge.eps \
            smooth-edge.eps \
            $(srcdir)/versenfuse.texi $(top_builddir)/config-h.texi \
            $(enfuse_TEXINFOS)

enfuse.html \
enfuse.pdf: photographic-workflow.png \
            focus-stack-decision-tree.png \
            entropy.png \
            entropy-cutoff.png \
            gaussian.png \
            laplacian-of-gaussian.png \
            local-analysis-window.png \
            sharp-edge.png \
            smooth-edge.png
