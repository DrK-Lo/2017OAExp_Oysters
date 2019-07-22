## Both packages should be installed ahead of time from CRAN
library(tximport)
library(readr)
library(DESeq2)

setwd("/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files")

### Sample Metadata (needed for DESeqData Object)
model<-read.csv("metadata/metadata_cvirginica_rna_meta.txt", header=TRUE)
model$Treatment <- as.factor(model$treatment)
model$Time <- as.character(model$timepoint)
model$Time[model$Time == "3"] <- "09"
model$Time[model$Time == "6"] <- "80"
model$Time <- as.factor(model$Time)
model$Pop <- as.factor(model$population)
model$Lane <- as.factor(model$lane)
model$SFV <-interaction(model$Time,model$Treatment) # Creates single factor variable for combi$

tranAggr <- readRDS("run20190610_counts/tranMatrixFull.RData")

print(head(tranAggr$counts))

dds_t <- DESeqDataSetFromTximport(tranAggr,
                                  model,
                                  design = ~ Lane + Pop + Treatment + Time + Treatment:Time)
print(str(dds_t))

saveRDS(dds_t,"tranDESeqDataObj.RData")

