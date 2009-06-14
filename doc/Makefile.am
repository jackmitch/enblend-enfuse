info_TEXINFOS = enblend.texi \
                enfuse.texi

enblend_TEXINFOS = auxmac.texi \
                   enblend-version.texi \
                   fdl.texi \
                   helpful-programs.texi \
                   tuning-memory-usage.texi \
                   understanding-masks.texi

enfuse_TEXINFOS = auxmac.texi \
                  enfuse-version.texi \
                  fdl.texi \
                  helpful-programs.texi \
                  tuning-memory-usage.texi \
                  understanding-masks.texi

#man_MANS = enblend.1 enfuse.1
#EXTRA_DIST = enblend.1 enfuse.1 enblend.texi enfuse.texi

# created by make(1), want to rebuild
MOSTLYCLEANFILES =

# created by make(1)
CLEANFILES = entropy.{txt,png,eps,pdf} \
             entropy-cutoff.{txt,png,eps,pdf} \
             gaussian.{txt,png,eps,pdf} \
             laplacian-of-gaussian.{txt,png,eps,pdf} \
             local-analysis-window.{txt,png,eps,pdf} \
             sharp-edge.{txt,png,eps,pdf} \
             smooth-edge.{txt,png,eps,pdf}

# created by configure(1)
DISTCLEANFILES = texinfo.tex

# e.g. *.info
MAINTAINERCLEANFILES =


#


%.txt %.png %.eps %.pdf: %.gp config.gp config-edge.gp
	$(GNUPLOT) $<


%.txt: %.fig
	$(SED) --expression='1,/---BEGIN-TEXT---/d' \
               --expression='/---END-TEXT---/,$$d' \
               --expression='s/^# \?//' < $< > $@


%.png: %.fig
	$(FIG2DEV) -L png $< $@


%.eps: %.fig
	$(FIG2DEV) -L eps $< $@


%.pdf: %.fig
	$(FIG2DEV) -L pdf $< $@

#


enfuse.info: entropy.txt \
             entropy-cutoff.txt \
             gaussian.txt \
             laplacian-of-gaussian.txt \
             local-analysis-window.txt \
             sharp-edge.txt \
             smooth-edge.txt
