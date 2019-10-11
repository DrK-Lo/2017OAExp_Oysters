I wrote some code to plot Alan's methylation data on chromosomes. 

* The code takes all the CpGs identified by Bismark and the average proportion methylated (beta) as calculated by Alan.
  * Note that a beta of 0 could be (a) 0% methylation or (b) no data.  Since the data was MBD enriched, it might be assumed 
that no data is 0% methylation.
  * The total counts matrix can be used to inform whether there was any coverage for that paritcular cpg

* The code then averages the beta across all individuals in all treatments, and then computes summary statistics for 100KB windows

Code can be found at markdown_files/compareGonadMantle and the [html output can be viewed here](https://htmlpreview.github.io/?https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/compareGonadMantle/CompareGonadMantle.html)

From Yaamini, we would like analogous data frames to these (search for these lines in the html to see what these files look like):

* mantle <- readRDS("~/Google Drive/LotterhosLab/2017-2020-Oyster OA Epigenetics/2017_CV_AdultExposureExp/Data/DNAm/Data/CG_unstranded_BetaMatrix.RData")
  * loci in columns and inidivudals in rows
  * proportion of reads that were called as methylated
* mantle_pos <- readRDS("~/Google Drive/LotterhosLab/2017-2020-Oyster OA Epigenetics/2017_CV_AdultExposureExp/Data/DNAm/Data/CG_unstranded_summaryTable.RData")
  * information about each locus
* mantle_counts <- readRDS("~/Google Drive/LotterhosLab/2017-2020-Oyster OA Epigenetics/2017_CV_AdultExposureExp/Data/DNAm/Data/CG_unstranded_TotalCountMatrix.RData")
  * total number of reads in that individual at that location
* mantle_info <- readRDS("~/Google Drive/LotterhosLab/2017-2020-Oyster OA Epigenetics/2017_CV_AdultExposureExp/Data/meta/metadata_20190811.Rdata")
  * information about each individual sample
  
  Alternatively, Yaamini could provide us with methylated and unmethylated C counts for each individual and locus, and we 
  could do the calculations
