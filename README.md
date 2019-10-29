## Working Repository for 2017 Oyster OA exposure and differential gene expression and methylation study

### Overview  
This repository was created to keep track of the code and preliminary analyses for the 2017 Oyster Acidification Experiment. RNA was extracted from 24 total samples, 12 samples from Timepoint 3 (June 13th) and 12 samples from Timepoint 6 (August 22/24th), with six samples from the two extreme treatments (400 and 2800 ppm) within each of these timepoints. 

### [Manuscript (in progress)](https://docs.google.com/document/d/1UTjTN_KC_exGVlf0I0UpntO7woBzGoJ3CKZel8WXobc/edit?ts=5bbf8c38)

### [Notebook](https://github.com/epigeneticstoocean/2017OAExp_Oysters/tree/master/notebook) 

## Analysis (in progress)

### [Transcriptomic Analysis](https://github.com/epigeneticstoocean/2017OAExp_Oysters/tree/master/markdown_files/Transcriptomic)
  * [01 - Creating Feature Count Matrices]()
  * [02 - Filtering and Normalization](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/Transcriptomic/02_AE17_RNA_featureFiltering_objectCreation.md)
  * [03 - Multivariate Analyses]()
  * [04 - Differential Expression]()
  * [05 - Biomineralization Genes](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/Transcriptomic/05_AE17_RNA_biomineralizationGenes.md)
  * [Extra Scripts](https://github.com/epigeneticstoocean/2017OAExp_Oysters/tree/master/markdown_files/Transcriptomic/RNAseq_Additional)
  * [Diagnostic Scripts](https://github.com/epigeneticstoocean/2017OAExp_Oysters/tree/master/markdown_files/Transcriptomic/Diagnostics)

### [DNA methylation Analysis](https://github.com/epigeneticstoocean/2017OAExp_Oysters/tree/master/markdown_files/DNAm)

### [Genomic References](https://github.com/epigeneticstoocean/2017OAExp_Oysters/tree/master/markdown_files/References)

### [Phenotypes](https://github.com/epigeneticstoocean/2017OAExp_Oysters/tree/master/markdown_files/Phenotype_Analysis)
  * [EPF Full Timeseries](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/Phenotype_Analysis/1_AE17_epf_Total.md)
  * [Calcification](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/Phenotype_Analysis/3A_AE17_calcification_FullAnalysis.md)
  * [Full Carbonate Chemistry - in supplemental](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/Phenotype_Analysis/2_AE17_epf_carbChem.md)

### [Water Chemistry](https://github.com/epigeneticstoocean/2017OAExp_Oysters/tree/master/markdown_files/WC) 

### Extra (outdated)

- [Read Count Break Down](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/extra/readAnalysis.md): Summarized raw read, alignment, and final gene count data for each individual. This file also contains a preliminary look at the unfiltered ```geneCount``` matrix. Specifically, the proportion of counts that fall in the upper 10,5,and1% of genes and how the counts of those genes compare between treatments.
- [Gene Count Matrix Comparison](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/extra/starReRun_17005_comparison.md): This was a brief comparisons between the geneCounts for sample 17005 generated during the original Original STAR Pipeline run (MAR-2018) and new STAR Pipeline run which used the same data and sample parameters.
- [Creating Oyster Transcriptome Reference Files](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/extra/transcriptomeReferenceFile_generation.md): Script for converting ```.fna``` files from the NCBI website into a dataframe that is stored as an ```.RData``` file that can be used for downstream steps. For example, this file is used with the ```txImport``` function to aggregate transcript quantification from the salmon mapper into genes.   
- [Salmon vs STAR Mapping Comparison](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/extra/STAR_Salmon_Mapping_Comparison.md): Looked at gene count matrix created by using the STAR mapping pipeline (created MAR-2018) vs. Salmon (MAY-2019). Found Salmon led to much higher counts (x100 times), that only correlated with the STAR mapping outputs by ~ 11%.  

### RNASeq Data Analysis Directory

The **figures\/** directory holds a figure(s) output during the initial filtering steps to visualizing the alignment of sequences relative to various parts of the genome (introns, exons, coding regions, etc.).

The **notebook\/** directory contains a log of all major data entries and a brief description of what was done.

The **input_files\/** directory hold various files that are used in the initial bioinformatics processing steps, such as config and adapter files, and the preliminary data analysis (metadata).

The **markdown_files\/** holds a markdown for the initial bioinformatics pipeline used to get from raw sequences to a count matrix and the R code used for preliminary analyses. The initial steps mostly followed the EecSeq pipeline with some small modifications.

The **results\/** directory holds the final counts matrix for all genes across all 24 individuals.

The **src\/** folder holds scripts used in initial filtering steps. The dDocent script was modified to suit the needs of our sequencing specifications.

