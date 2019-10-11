
setwd("")

source("DNAmRefCode.R")
gene <- readRDS("")
exon <-  readRDS("")
meta <- readRDS("")

eMatch <- matchRef(meta,exon)
gMatch <- matchRef(meta,gene)
full_meta <-  cbind(meta,gMatch[,c(3,4,5)],eMatch[,c(4,5,9,10,11,12)])
colnames(full_meta) <- c("chr","pos","strand","motif","tri","gene","g_start","g_end","e_start","e_end","exon_ID","Parent","gene_ID","transcript_ID")

saveRDS(full_meta,"")