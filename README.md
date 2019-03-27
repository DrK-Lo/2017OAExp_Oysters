# Manuscript
[Manuscript](https://docs.google.com/document/d/1UTjTN_KC_exGVlf0I0UpntO7woBzGoJ3CKZel8WXobc/edit?ts=5bbf8c38)

## RNASeq Data Analysis Directory

This repository was created to keep track of the code and preliminary analyses for the 2017 Oyster Acidification Experiment. RNA was extracted from 24 total samples, 12 samples from Timepoint 3 (June 13th) and 12 samples from Timepoint 6 (August 22/24th), with six samples from the two extreme treatments (400 and 2800 ppm) within each of these timepoints. 

The **figures\/** directory holds a figure(s) output during the initial filtering steps to visualizing the alignment of sequences relative to various parts of the genome (introns, exons, coding regions, etc.).

The **input_files\/** directory hold various files that are used in the initial bioinformatics processing steps, such as config and adapter files, and the preliminary data analysis (metadata).

The **markdown_files\/** holds a markdown for the initial bioinformatics pipeline used to get from raw sequences to a count matrix and the R code used for preliminary analyses. The initial steps mostly followed the EecSeq pipeline with some small modifications.

The **results\/** directory holds the final counts matrix for all genes across all 24 individuals.

The **src\/** folder holds scripts used in initial filtering steps. The dDocent script was modified to suit the needs of our sequencing specifications.

