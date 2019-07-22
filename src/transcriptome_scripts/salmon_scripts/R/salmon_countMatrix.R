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
         countsFromAbundance = "dtuScaledTPM",
         txOut = TRUE,
         tx2gene = tr2gene)

print("Complete transcript counts...")

saveRDS(tranAggr,"tranMatrixFull.RData",compress = TRUE)
saveRDS(tranAggr$abundance,"tranMatrixAbundance.RData",compress = TRUE)
saveRDS(tranAggr$countsFromAbundance,"tranMatrixCountsFromAbundance.RData",compress = TRUE)
saveRDS(tranAggr$length,"tranMatrixLength.RData",compress = TRUE)

print("Creating DESeqDataSet....")

dds_g <- DESeqDataSetFromTximport(geneAggr,
                                  model,
                                  design = ~ Lane + Pop + Treatment + Time + Treatment:Time)

dds_t <- DESeqDataSetFromTximport(tranAggr,
                                  model,
                                  design = ~ Lane + Pop + Treatment + Time + Treatment:Time)

saveRDS(dds_g,"geneDESeqDataObj.RData")
saveRDS(dd_t,"tranDESeqDataObj.RData")

print("Successfully saved outputs...")
