# Performing differential methylation analysis using a bayesian mixed effects binomial model approach 

**Description**
After conversation with Katie we decided to go with a bayesian binomial mixed effects modelling approach implemented using the packaged BRMS. This operates as a wrapper around an independent coding language called `stan`, and provides `lme4` based syntax using bayesian model approach. This was selected because we were experiencing issues with model convergence when using a non-bayesian approach implemented in `GLMM` using a similarly structured binomial mixed effects model. 

**Input Data Notes**
I wrote code to perform this model on a provided matrix of beta values (methylC/totalC) for CpGs with at least 5x coverage per sample and located in genic regions (N=294261). Importantly, these counts are based on stranded data. so only consider coverage information for the cytosine located on the strand with the gene.

**Performing the model**
The model takes a long time to run so I wrote two scripts (UNIX and R) to peform the regression and store the results. The [R script](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/DNAm/cluster_scripts/R/07_diffMethylation_brms.R) does the heavy lifting, calling the BRMS functions, fitting the model, and saving the outputs. This should require minimal (if any change) between runs. Next, I have a [shell script](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/DNAm/cluster_scripts/07_diffMethylation_brms.sh) that allows for the pseudo parralelization of the R script, by dividing the specified range of loci being analyzed into even blocks for the number of cores you specify. It will then call the R script that many times and run it in the background until it completes. 

Things you may need to change in shell script:
* `DIR="/shared_lab/20180226_RNAseq_2017OAExp/DNAm"`
* `INPUT="processed_samples/05_countSummary"` - path from directory to input data folder
* `RSCRIPT="scripts/R/"` - path from your directory to the R script
* `OUTPUT="processed_samples/07_brmsSummary"` - path from directory to output folder
* `DIAG="/diagnostic/"` - nested folder in output were diagnostic results go 
* `START=1` - Starting cpg position (row) within the input data frame
* `END=100000` - Ending cpg position (row) within the input data frame
* `NCORES=50` - Number of cores 

*Note* - Code not error proof - the difference between START and END should be greater than then number of cores

**Output Data Notes**
* `/DNAm_gene_BRMS_",MINp,"_",MAXp,"_modelParam.RData` 
  * List or model parameters -  these will be the same for each shell script run.
* `/DNAm_gene_BRMS_",MINp,"_",MAXp,"_modelSummary.RData`
  * Model summary stats as list
* `/DNAm_gene_BRMS_",MINp,"_",MAXp,"_modelMarginalEffects.RData`
  * Marginal effects as list
* `/DNAm_gene_BRMS_",MINp,"_",MAXp,"_plannedComparisons.RData`
  * Primary and planned comparison outputs as list
  
