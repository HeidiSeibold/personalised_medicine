ORGDIR = ~/Documents/svn/WhatIf/WhatIf_forests/analysis

BASEFILES = \
	basis/personalised_models.R \
	basis/variable_importance.R \
	basis/dependence_plots.R \
	$(DONE)

INFOS = \
	personalised_models.Rnw \
	log_likelihoods.Rnw \
	variable_importance.Rnw \
	partial_dependence_plots.Rnw \
	$(DONE)

RNWFILES = \
	   simulation.Rnw \
	   ALSFRS.Rnw \
	   ALSsurv.Rnw \
	   $(DONE)

OTHER = \
	ALSFRS_R.R \
	$(INFOS) \
	ALSsurv_R.R \
	$(DONE)


# copying files from svn project
basis: $(ORGDIR)/$(BASEFILES)
	rsync -rt $(ORGDIR)/basis/ ./basis/ 

$(RNWFILES): $(ORGDIR)/$@
	rsync -t $(ORGDIR)/$@ ./$@

$(OTHER): $(ORGDIR)/$@
	rsync -t $(ORGDIR)/$@ $@
	
#.PHONY: copy
copy: $(RNWFILES) basis $(OTHER)


# run Rnw files
simulation.pdf: $(BASEFILES) $(INFOS) simulation.Rnw
	R CMD Sweave simulation.Rnw
	pdflatex simulation.tex

ALSFRS.pdf: $(BASEFILES) ALSFRS.Rnw ALSFRS_R.R  
	R CMD Sweave ALSFRS.Rnw
	pdflatex ALSFRS.tex
	R CMD Sweave ALSFRS_figures.Rnw
	pdflatex ALSFRS_figures.tex

ALSsurv.pdf: $(BASEFILES) ALSsurv.Rnw ALSsurv_R.R
	R CMD Sweave ALSsurv.Rnw
	pdflatex ALSsurv.tex
	R CMD Sweave ALSsurv_figures.Rnw
	pdflatex ALSsurv_figures.tex

clean:
	-rm -rf *.RData
	-rm -f *.out
	-rm -f *.log
	-rm -f *.aux
	-rm -f *.bbl
	-rm -f *.blg
	-rm -f ALSFRS.tex
	-rm -f simulation.tex
	-rm -f ALSsurv.tex
	-rm -f ALSFRS.R
	-rm -f ALSsurv.R
	-rm -f simulation.R


