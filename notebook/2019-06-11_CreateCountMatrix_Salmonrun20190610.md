# Create count matrices from the salmon `quant.sf` outputs from the `run20190610` salmon run

**Description**: To convert the transcript quantity outputs from salmon into either transcript or gene (isoform) level counts I used the R package `tximport`, which will aggregate (or not in the case of the the transcripts) counts of a set of quant.sf files (one for each sample) and arrange them into a convienent matrix format (Samples x Loci). Furthermore, I also used this package to scale my counts by both the length of the isoform and the size of the library (using `countsFromAbundance = "lengthScaledTPM"` arguement). These matrices can then be used for either DE analysis or other downstream analysis.

[LINK to `tximport` website](http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)

**Input Info**

`quant.sf` files for each sample in the `run20190610` run
* PATH : `/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files/run20190610/`

Transcript list (based on converted .fna transcriptome)
* PATH : `/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files/run20190610/NCBI_transcriptome/transcriptome_table.RData`

Metadata (used for created DESeqData Object)
* PATH : `/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files/run20190610/metadata/metadata_cvirginica_rna_meta.txt'

**Output** 

* All outputs were saved at `/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files/run20190610_counts`

Output 1: The default object created by `tximport`
* Files
  * Transcripts: `tranMatrixFull.RData`
  * Genes: `geneMatrixFull.RData`
* These files contains multiple pieces of data including:
  * A matrix of 'raw' counts : `$counts`
  * A matrix of adbundance estimates : `$abundance`
  * A vector of transcript lengths : `$length`
  * A matrix of counts scaled : `$countsFromAbundance`

Output 2: Feature Abundance Matrix ($abundance)
* Files
  * Transcripts: `tranMatrixAbundance.RData`
  * Genes: `geneMatrixAbundance.RData`

Output 3: Feature Counts from Abundance Matrix ($countsFromAbundance)
* Files 
  * Transcripts: 'tranMatrixCountsFromAbundance.Rdata'
  * Genes: 'geneMatrixCountsFromAbundance.Rdata'

Output 4: Transcript Lengths ($length)
* Files 
  * Transcripts: 'tranMatrixCountsLength.Rdata'
  * Genes: 'geneMatrixCountsLength.Rdata'

Output 5: DESeqData Object (for use in DESeq2 DE analysis)
* Files
  * Transcripts: 'tranDESeqDataObj.RData'
  * Genes : 'geneDESeqDataObj.RData'
  
**Additional Thoughts**
* Output 1 (the full output) is quite large and does not want to open on my local machine (can't allocate enough memory), this prompted me to save the different matrices separately (Abundance,CountsFromAbundance,Length).

* I also created the DESeq Data object here to see if this was a smaller size and might be more management for local computing resources. (I think one reason the full matrix is so large is that it saved the inferential replicates).

* The authors of `tximport` strongly recommend only using the counts from abundance estimates for downstream analysis. The rationale being: "Simply passing the summed estimated transcript counts does not correct for potential differential isoform usage (the offset), which is the point of the tximport methods (Soneson, Love, and Robinson 2015) for gene-level analysis" (from tximport website)

* I also ran generated inferential replicates for each transcript during my salmon run, these are incorporated into the tximport function and stored in the full object. I don't think this impacts the count estimates, but it is possible to examine gene level variance within the created object.

### Code

**Command line**
```
Rscript 01_salmon_countMatrix.R
```
NOTE: you have to be within `/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files/scripts` folder for this to work, or set the full path within you call the `.R` script.

**R script for `01_salmon_countMatrix.R`**:

```
## This is the general remote cluster version of the 01_salmon_countmatrix step.
## See https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/Salmon_pipeline/README.md
## for additional details.

## Both packages should be installed ahead of time from CRAN
library(tximport)
library(readr)
library(DESeq2)

setwd("/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files")
## NOTE: be sure the both the folder with your salmon transcript quantifications and the file with the transcript summaries are located within this path

### Directory of sample salmon files
DIR <- "run20190610/" # SET THIS to the sample directory with the files you want included in gene count matrix
ls_files <- list.files(DIR)
files_input <- file.path(DIR,ls_files,"/quant.sf",fsep = "")

print("Files being run...")
print(files_input)


### Transcript List
TRANS <- readRDS("NCBI_transcriptome/transcriptome_table.RData") # SET THIS to name of file with your list of transcripts
tr2gene <- data.frame(TXNAME = TRANS$fullID,GENEID=TRANS$location)
tr2gene$TXNAME <- as.character(tr2gene$TXNAME)
tr2gene$GENEID <- as.character(tr2gene$GENEID)

### Sample Metadata (needed for DESeqData Object)
model<-read.csv("metadata/metadata_cvirginica_rna_meta.txt", header=TRUE)
model$Treatment <- as.factor(model$treatment)
model$Time <- as.character(model$timepoint)
model$Time[model$Time == "3"] <- "09"
model$Time[model$Time == "6"] <- "80"
model$Time <- as.factor(model$Time)
model$Pop <- as.factor(model$population)
model$Lane <- as.factor(model$lane)
model$SFV <-interaction(model$Time,model$Treatment) # Creates single factor variable for combination of time and treatment

print("Directory set and transcript list configure...")

### RUN TXIMPORT

print("Starting counting...")

# For gene level quantification
geneAggr <- tximport(files_input,
            type = "salmon",
            countsFromAbundance = "lengthScaledTPM",
            tx2gene = tr2gene)

# Save gene files
saveRDS(geneAggr,"geneMatrixFull.RData",compress = TRUE)
saveRDS(geneAggr$abundance,"geneMatrixAbundance.RData",compress = TRUE)
saveRDS(geneAggr$countsFromAbundance,"geneMatrixCountsFromAbundance.RData",compress = TRUE)
saveRDS(geneAggr$length,"geneMatrixLength.RData",compress = TRUE)

print("Completed gene counts...")


# For transcript level quantification
tranAggr <- tximport(files_input,
         type = "salmon",
         countsFromAbundance = "lengthScaledTPM",
         txOut = TRUE)

print("Complete transcript counts...")

saveRDS(tranAggr,"tranMatrixFull.RData",compress = TRUE)
saveRDS(tranAggr$abundance,"tranMatrixAbundance.RData",compress = TRUE)
saveRDS(tranAggr$countsFromAbundance,"tranMatrixCountsFromAbundance.RData",compress = TRUE)
saveRDS(tranAggr$length,"tranMatrixLength.RData",compress = TRUE)

print("Creating DESeqDataSet....")

dds_g <- DESeqDataSetFromTximport(geneAggr,
                                  model,
                                  design = ~ Lane + Pop + Treatment + Time + Treatment:Time)

dds_t <- DESeqDataSetFromTximport(geneAggr,
                                  model,
                                  design = ~ Lane + Pop + Treatment + Time + Treatment:Time)

saveRDS(dds_g,"geneDESeqDataObj.RData")
saveRDS(dd_t,tranDESeqDataObj.RData")

print("Successfully saved outputs...")
```
