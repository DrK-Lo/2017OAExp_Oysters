## Working Repository for 2017 Oyster OA exposure and differential gene expression and methylation study

### Overview  
This repository was created to keep track of the code and preliminary analyses for the 2017 Oyster Acidification Experiment. RNA was extracted from 24 total samples, 12 samples from Timepoint 3 (June 13th) and 12 samples from Timepoint 6 (August 22/24th), with six samples from the two extreme treatments (400 and 2800 ppm) within each of these timepoints. 

### [Manuscript (in progress)](https://docs.google.com/document/d/1UTjTN_KC_exGVlf0I0UpntO7woBzGoJ3CKZel8WXobc/edit?ts=5bbf8c38)

### [Notebook](https://github.com/epigeneticstoocean/2017OAExp_Oysters/tree/master/notebook) 

### Main Analysis (in progress)

STAR Pipeline
  - [Gene Matrix Filtering](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/STAR_pipeline/03A_CV17_RNA_countFilteringandAnalysis.md)
  - [Gene Normalization](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/STAR_pipeline/03B_CV17_RNA_countAnalysis.md)
  - [DAPC Analysis](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/STAR_pipeline/04B_CV17_RNA_DAPC.md)
  - [Targeted Gene Queries](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/STAR_pipeline/04C_CV17_RNA_targetGeneQuery.md)
  - [RDA, DAPC, PERMANOVA(adonis), other Multivariate visualize and statistics](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/STAR_pipeline/04D_CV17_RNA_CCAandRDA.md)
  - [DESeq Analysis for Differential Expression](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/STAR_pipeline/04E_CV17_RNA_DESeqAnalysis.md)

[SALMON Pipeline](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/Salmon_pipeline/README.md)
  - [Transcript Quantification]()
  - [RDA, nMDS, PERMANOVA(adonis), other Multivariate visualize and statistics]()
  - [Differential Expression - DESeq2]()
  - [Differential Expression - Sleuth Analysis]()
  - [Correlated Gene Networks]()

Phenotypes 
  - [EPF Phenotype Analysis](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/Phenotype_Analysis/AE17_epfPhenotype.md)

[Water Chemistry]() 


### Extra 

- [Read Count Break Down](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/extra/readAnalysis.md): Summarized raw read, alignment, and final gene count data for each individual. This file also contains a preliminary look at the unfiltered ```geneCount``` matrix. Specifically, the proportion of counts that fall in the upper 10,5,and1% of genes and how the counts of those genes compare between treatments.
- [Gene Count Matrix Comparison](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/extra/starReRun_17005_comparison.md): This was a brief comparisons between the geneCounts for sample 17005 generated during the original Original STAR Pipeline run (MAR-2018) and new STAR Pipeline run which used the same data and sample parameters.
- [Salmon vs STAR Mapping Comparison](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/extra/STAR_Salmon_Mapping_Comparison.md): Looked at gene count matrix created by using the STAR mapping pipeline (created MAR-2018) vs. Salmon (MAY-2019). Found Salmon lead to much higher counts (x100 times), that only correlated with the STAR mapping outputs by ~ 11%.  

### RNASeq Data Analysis Directory

The **figures\/** directory holds a figure(s) output during the initial filtering steps to visualizing the alignment of sequences relative to various parts of the genome (introns, exons, coding regions, etc.).

The **notebook\/** directory contains a log of all major data entries and a brief description of what was done.

The **input_files\/** directory hold various files that are used in the initial bioinformatics processing steps, such as config and adapter files, and the preliminary data analysis (metadata).

The **markdown_files\/** holds a markdown for the initial bioinformatics pipeline used to get from raw sequences to a count matrix and the R code used for preliminary analyses. The initial steps mostly followed the EecSeq pipeline with some small modifications.

The **results\/** directory holds the final counts matrix for all genes across all 24 individuals.

The **src\/** folder holds scripts used in initial filtering steps. The dDocent script was modified to suit the needs of our sequencing specifications.

