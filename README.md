# Personalised Medicine

This repository contains the code to our paper on personalised medicine.

We propose model-based random forests as a method to detect
similarities between patients with respect to their treatment effect
and on this basis compute personalised models for new patients to
predict their individual treatment effect.  We apply this method to a
database of several clinical trials on ALS and show that some patients
have a systematically higher benefit from the drug Riluzole than
others. 

We use the [partykit](https://cran.r-project.org/web/packages/partykit) package (especially functions `cforest()` and `predict.cforest()`) and parametric 
models such as `survreg()` from package [survival](https://cran.r-project.org/web/packages/survival) or `glm()`.
In order to run the code please download the data from https://nctu.partners.org/ProACT and run the
corresponding demo from the R package [TH.data](https://cran.r-project.org/web/packages/TH.data). Store the data in a subfolder called /data.


## Background and aim
Medical breakthroughs have become rare due to that fact that the problems left
to solve usually are highly complex and diseases may have unknown sub-diseases
that demand different treatments. Clinical trials fail at finding treatments
that are effective for a broad group of patients. Thus it becomes more and more
important to make individualised treatment decisions.

Amyotrophic lateral sclerosis (ALS) is a complex disease and the only approved
drug to treat the disease is Riluzole. Riluzole, however, only has a minor
impact on prolonging the survival of the average patient 
([European Medicines Agency, 2012](http://www.ema.europa.eu/docs/en_GB/document_library/EPAR_-_Summary_for_the_public/human/002622/WC500127609.pdf)) 
and often comes along with side effects.  Some ALS patients might
have a greater benefit from Riluzole than suggested by the overall treatment
effect and it is important for a new patient to know if she or he is likely to
have a reasonable effect or if the treatment effect is not strong enough to
tolerate the common side effects.

We developed the proposed method not only to identify patient characteristics
that should be considered in treatment decisions but also and primarily to
predict the treatment effect of a future patient and thus enhance treatment
decisions.

We apply our method to the PRO-ACT database 
([Atassi et al, 2014](http://dx.doi.org/10.1212/WNL.0000000000000951)) containing
fully de-identified data of ALS patients from several clinical trials. We 
are primarily interested in the effect of Riluzole on the survival of ALS
patients. 


## Methods
Model-based recursive partitioning can be used as a method to detect subgroups
with differing treatment effects (Seibold et al, 2015 :clock5:).  The tree is grown
based on the correlation between patient characteristics and the score function
of a given base model (see [Zeileis et al, 2008](http://dx.doi.org/10.1198/106186008X319331)).  
In the context of subgroup
detection for differential treatment effects the primary endpoint in the base
model is some health measurement or the survival time and the covariate is the
(usually binary) treatment indicator.  In this application the primary outcome
is the survival time of patients in days from treatment start and as base model
we use a Weibull model.  

Rather than giving treatment effects within subgroups - as is the case in
model-based recursive partitioning - we want go one step further and estimate
personalised treatment effects here. To be able to do this, in a first step we
compute random forests based on multiple model-based trees t=1,...,T. Each
tree is based on a subsample of the learning data and a
subsample of eligible patient characteristics that possibly define the
subgroups.  The forest is used merely as a method to detect similarity between
patients with respect to expected outcome and, more importantly, treatment
effect. So in a second step we compute the similarity between patients.  The
similarity of a patient i to any patient given in the learning data is
defined as the number of times the patients are assigned to the same subgroup
in the trees ([Hothorn et al, 2004](http://dx.doi.org/10.1002/sim.1593)).  
Finally the expected treatment effect for this
new patient is computed using the base model and a weighted version of the
learning data.  The weights that are used for the model are the similarities.
In this way one can compute a separate model and thus a separate treatment
effect for every patient. 
Dependence plots, i.e., plotting the expected outcome improvement by the
treatment against the single patient characteristics, give us an idea on the
connection between the two. 




:clock5: *Seibold, H., Zeileis, A. and Hothorn, T. (accepted 2015-10-02). Model-based
recursive partitioning for subgroup analyses. International Journal of Biostatistics.*
