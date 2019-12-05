README

For visualizing DNA methylation interactions, significant CpGs were clustered into groups based on their scaled treatment means. Treatment means were scaled relative to the overall mean across 4 treatments for that CpG, so that whether the treatment tended to increase or decrease methylation could be visualized.

# DATA

The "mainEffects.." file is a matrix of ones and zeros, where 1 means that the CpG was significant for either time, 
treatment, or the interaction (depending on the column). 
The first column is an index that refers to the specific row in the "MeanMethylation" file (which is sorted). 

So if you want all loci that were significant for any effect you could simply use all the indexes from "mainEffects.." 
to subset the "MeanMethylation..." file.

CpGMeanByTreatment_withSignificantInteraction - only ones with significant interactions

# FIGS

#### Sig Codes:

* 1 = Main effect of time, no interaction
* 2 = Main effect of treatment, no interaction
* 3 = Main effect of time and treatment, no interaction (no loci)
* 4 = Interaction between treatment and time

#### Chrom plots:
Chromosome location was scaled between 0 and 1 for each chromosome


