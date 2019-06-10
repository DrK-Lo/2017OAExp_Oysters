
## This is the general remote cluster version of the 01_salmon_countmatrix step.
## See https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/Salmon_pipeline/README.md
## for additional details.

## Both packages should be installed ahead of time from CRAN
library(tximport)
library(readr)

setwd("/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files")
## NOTE: be sure that both the folder with your salmon transcript quantifications and the file with the transcript summaries are located within this path

### Directory of sample salmon files
DIR <- "/run20180512/" # SET THIS to the sample directory with the files you want included in gene count matrix
ls_files <- list.files(DIR)
files_input <- file.path(DIR,ls_files,"/quant.sf",fsep = "")

### Transcript List 
TRANS <- readRDS("transcriptome_fromGenome_table.RData") # SET THIS to name of file with your list of transcripts
tr2gene <- data.frame(TXNAME = TRANS$fullID,GENEID=TRANS$location)
tr2gene$TXNAME <- as.character(tr2gene$TXNAME)
tr2gene$GENEID <- as.character(tr2gene$GENEID)


### RUN TXIMPORT  
  
# For gene aggregation  
geneAggr <- tximport(files_input,
            type = "salmon",
            tx2gene = tr2gene)
saveRDS(geneAggr,"/run20180512_matrixOutputs/geneMatrixFull_default.RData",compress = TRUE)
saveRDS(geneAggr$abundance,"/run20180512_matrixOutputs/geneMatrixAbundance_default.RData",compress = TRUE)

# For transcript level quantification
#NOTE additional arguement called countsFromAbundance
# Scales counts using the average transcript length over the samples
# and then the library size.

tranAggr <- tximport(files_input,
         type = "salmon",
         countsFromAbundance = "lengthScaledTPM",
         txOut = TRUE)
saveRDS(tranAggr,"/run20180512_matrixOutputs/tranMatrixFull_default.RData",compress = TRUE)
saveRDS(tranAggr$abundance,"/run20180512_matrixOutputs/tranMatrixAbundance_default.RData",compress = TRUE)


