ORGDIR = ~/Documents/svn/WhatIf/WhatIf_forests/analysis/

BASEFILES = \
	basis/personalised_models.R \
	basis/pruning.R \
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
	   simulation2.Rnw \
	   ALSFRS.Rnw \
	   ALSsurv.Rnw \
	   $(DONE)

OTHER = \
	ALSFRS_R.R \
	$(INFOS) \
	ALSsurv_R.R \
	$(DONE)



# copying files from svn project
basis: $(addprefix $(ORGDIR),$(BASEFILES))
	rsync -rt $(addprefix $(ORGDIR),$(BASEFILES)) ./basis/ 

$(RNWFILES): $(addprefix $(ORGDIR),$(RNWFILES))
	rsync -t  $(addprefix $(ORGDIR),$(RNWFILES)) .

$(OTHER): $(addprefix $(ORGDIR),$(OTHER))
	rsync -t $(addprefix $(ORGDIR),$(OTHER)) .

#.PHONY: test
test: 
	echo $(addprefix $(ORGDIR),$(RNWFILES))
	
#.PHONY: copy
copy: $(RNWFILES) basis $(OTHER)


# run Rnw files
simulation.pdf: $(BASEFILES) $(INFOS) simulation.Rnw
	R CMD Sweave simulation.Rnw
	pdflatex simulation.tex

simulation2.pdf: $(BASEFILES) $(INFOS) simulation2.Rnw
	R CMD Sweave simulation2.Rnw
	pdflatex simulation2.tex


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


