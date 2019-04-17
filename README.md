## Working Repository for 2017 Oyster OA exposure and differential gene expression and methylation study

### Overview  
This repository was created to keep track of the code and preliminary analyses for the 2017 Oyster Acidification Experiment. RNA was extracted from 24 total samples, 12 samples from Timepoint 3 (June 13th) and 12 samples from Timepoint 6 (August 22/24th), with six samples from the two extreme treatments (400 and 2800 ppm) within each of these timepoints. 

### Manuscript (in progress)
[Manuscript](https://docs.google.com/document/d/1UTjTN_KC_exGVlf0I0UpntO7woBzGoJ3CKZel8WXobc/edit?ts=5bbf8c38)

### Main Analysis (in progress)

- [Gene Matrix Filtering](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/STAR_pipeline/03A_CV17_RNA_countFilteringandAnalysis.md)
- [EPF Phenotype Analysis - In progress](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/Phenotype_Analysis/AE17_epfPhenotype.md)

### Extra 

- [Read Count Break Down](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/extra/readAnalysis.md): Summarized raw read, alignment, and final gene count data for each individual. This file also contains a preliminary look at the unfiltered ```geneCount``` matrix. Specifically, the proportion of counts that fall in the upper 10,5,and1% of genes and how the counts of those genes compare between treatments.

### RNASeq Data Analysis Directory

The **figures\/** directory holds a figure(s) output during the initial filtering steps to visualizing the alignment of sequences relative to various parts of the genome (introns, exons, coding regions, etc.).

The **notebook\/** directory contains a log of all major data entries and a brief description of what was done.

The **input_files\/** directory hold various files that are used in the initial bioinformatics processing steps, such as config and adapter files, and the preliminary data analysis (metadata).

The **markdown_files\/** holds a markdown for the initial bioinformatics pipeline used to get from raw sequences to a count matrix and the R code used for preliminary analyses. The initial steps mostly followed the EecSeq pipeline with some small modifications.

The **results\/** directory holds the final counts matrix for all genes across all 24 individuals.

The **src\/** folder holds scripts used in initial filtering steps. The dDocent script was modified to suit the needs of our sequencing specifications.

