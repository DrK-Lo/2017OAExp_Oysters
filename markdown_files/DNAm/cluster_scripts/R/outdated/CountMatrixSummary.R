#### Libraries ####
library(dplyr)
library(ggplot2)
library(matrixStats)
#library(R.utils)

sprintf("Reading in data")
#### Read in Data ####

setwd("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/processed_samples")
# Input / Output Folders
INPUT <- ("04_countMatrix/")
OUTPUT <- ("05_countSummary/")

gene <- readRDS("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/reference/gene_GeneLoc.RData")
exon <- readRDS("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/reference/Exon_GeneLoc.RData")

## Methylation Data
# Stranded
sumTab_stranded <- readRDS(paste0(INPUT,"CG_stranded_summaryTable.RData"))
umC_stranded <-    readRDS(paste0(INPUT,"CG_stranded_unmethylCountMatrix.RData"))
mC_stranded <-     readRDS(paste0(INPUT,"CG_stranded_methylCountMatrix.RData"))
tC_stranded <-     readRDS(paste0(INPUT,"CG_stranded_TotalCountMatrix.RData"))
beta_stranded <-   readRDS(paste0(INPUT,"CG_stranded_BetaMatrix.RData"))
# Unstranded
sumTab_unstranded <- readRDS(paste0(INPUT,"CG_unstranded_summaryTable.RData"))
umC_unstranded <-    readRDS(paste0(INPUT,"CG_unstranded_unmethylCountMatrix.RData"))
mC_unstranded <-     readRDS(paste0(INPUT,"CG_unstranded_methylCountMatrix.RData"))
tC_unstranded <-     readRDS(paste0(INPUT,"CG_unstranded_TotalCountMatrix.RData"))
beta_unstranded <-   readRDS(paste0(INPUT,"CG_unstranded_BetaMatrix.RData"))

sprintf("Pruning by min count threshold...")
# All cytosine
# Pruning by variable min count thresholds
tC_min_u <- rowMins(tC_unstranded) # vector of minimum counts for each cpg
tC_min_s <- rowMins(tC_stranded) # vector of minimum counts for each cpg
### At least 1 count per sample
#At least 1 count per sample (unstranded)
meta_u_1 <- sumTab_unstranded[tC_min_u>=1,]
mC_u_1 <- mC_unstranded[tC_min_u>=1,]
umC_u_1 <- umC_unstranded[tC_min_u>=1,]
tC_u_1 <- tC_unstranded[tC_min_u>=1,]
beta_u_1 <- beta_unstranded[tC_min_u>=1,]
#At least 1 count per sample (stranded)
meta_s_1 <- sumTab_stranded[tC_min_s>=1,]
mC_s_1 <- mC_stranded[tC_min_s>=1,]
umC_s_1 <- umC_stranded[tC_min_s>=1,]
tC_s_1 <- tC_stranded[tC_min_s>=1,]
beta_s_1 <- beta_stranded[tC_min_s>=1,]
#Saving unstranded
saveRDS(meta_u_1,paste0(OUTPUT,"CG_unstranded_summaryTable_1.RData"))
saveRDS(mC_u_1,paste0(OUTPUT,"CG_unstranded_mC_1.RData"))
saveRDS(umC_u_1,paste0(OUTPUT,"CG_unstranded_umC_1.RData"))
saveRDS(tC_u_1,paste0(OUTPUT,"CG_unstranded_tC_1.RData"))
saveRDS(beta_u_1,paste0(OUTPUT,"CG_unstranded_beta_1.RData"))
#Saving stranded
saveRDS(meta_s_1,paste0(OUTPUT,"CG_stranded_summaryTable_1.RData"))
saveRDS(mC_s_1,paste0(OUTPUT,"CG_stranded_mC_1.RData"))
saveRDS(umC_s_1,paste0(OUTPUT,"CG_stranded_umC_1.RData"))
saveRDS(tC_s_1,paste0(OUTPUT,"CG_stranded_tC_1.RData"))
saveRDS(beta_s_1,paste0(OUTPUT,"CG_stranded_beta_1.RData"))

### At least 5 count per sample
#At least 5 count per sample (unstranded)
meta_u_5 <- sumTab_unstranded[tC_min_u>=5,]
mC_u_5 <- mC_unstranded[tC_min_u>=5,]
umC_u_5 <- umC_unstranded[tC_min_u>=5,]
tC_u_5 <- tC_unstranded[tC_min_u>=5,]
beta_u_5 <- beta_unstranded[tC_min_u>=5,]
#At least 5 count per sample (stranded)
meta_s_5 <- sumTab_stranded[tC_min_s>=5,]
mC_s_5 <- mC_stranded[tC_min_s>=5,]
umC_s_5 <- umC_stranded[tC_min_s>=5,]
tC_s_5 <- tC_stranded[tC_min_s>=5,]
beta_s_5 <- beta_stranded[tC_min_s>=5,]
#Saving unstranded
saveRDS(meta_u_5,paste0(OUTPUT,"CG_unstranded_summaryTable_5.RData"))
saveRDS(mC_u_5,paste0(OUTPUT,"CG_unstranded_mC_5.RData"))
saveRDS(umC_u_5,paste0(OUTPUT,"CG_unstranded_umC_5.RData"))
saveRDS(tC_u_5,paste0(OUTPUT,"CG_unstranded_tC_5.RData"))
saveRDS(beta_u_5,paste0(OUTPUT,"CG_unstranded_beta_5.RData"))
#Saving stranded
saveRDS(meta_s_5,paste0(OUTPUT,"CG_stranded_summaryTable_5.RData"))
saveRDS(mC_s_5,paste0(OUTPUT,"CG_stranded_mC_5.RData"))
saveRDS(umC_s_5,paste0(OUTPUT,"CG_stranded_umC_5.RData"))
saveRDS(tC_s_5,paste0(OUTPUT,"CG_stranded_tC_5.RData"))
saveRDS(beta_s_5,paste0(OUTPUT,"CG_stranded_beta_5.RData"))

### At least 10 count per sample
#At least 10 count per sample (unstranded)
meta_u_10 <- sumTab_unstranded[tC_min_u>=10,]
mC_u_10 <- mC_unstranded[tC_min_u>=10,]
umC_u_10 <- umC_unstranded[tC_min_u>=10,]
tC_u_10 <- tC_unstranded[tC_min_u>=10,]
beta_u_10 <- beta_unstranded[tC_min_u>=10,]
#At least 10 count per sample (stranded)
meta_s_10 <- sumTab_stranded[tC_min_s>=10,]
mC_s_10 <- mC_stranded[tC_min_s>=10,]
umC_s_10 <- umC_stranded[tC_min_s>=10,]
tC_s_10 <- tC_stranded[tC_min_s>=10,]
beta_s_10 <- beta_stranded[tC_min_s>=10,]
#Saving unstranded
saveRDS(meta_u_10,paste0(OUTPUT,"CG_unstranded_summaryTable_10.RData"))
saveRDS(mC_u_10,paste0(OUTPUT,"CG_unstranded_mC_10.RData"))
saveRDS(umC_u_10,paste0(OUTPUT,"CG_unstranded_umC_10.RData"))
saveRDS(tC_u_10,paste0(OUTPUT,"CG_unstranded_tC_10.RData"))
saveRDS(beta_u_10,paste0(OUTPUT,"CG_unstranded_beta_10.RData"))
#Saving stranded
saveRDS(meta_s_10,paste0(OUTPUT,"CG_stranded_summaryTable_10.RData"))
saveRDS(mC_s_10,paste0(OUTPUT,"CG_stranded_mC_10.RData"))
saveRDS(umC_s_10,paste0(OUTPUT,"CG_stranded_umC_10.RData"))
saveRDS(tC_s_10,paste0(OUTPUT,"CG_stranded_tC_10.RData"))
saveRDS(beta_s_10,paste0(OUTPUT,"CG_stranded_beta_10.RData"))


# Creates a single geneToExonIndex file
# geneToExonIndex <- match(exon$gene_id,gene$gene_id)
# sum(is.na(geneToExonIndex)) # 0 confirms there is matching gene for each exon
# gene_M <- gene[geneToExonIndex,]
# gene_rename <- c("gene_chr","gene_method","gene_feature",
#                  "gene_start","gene_end","gene_score","gene_strand",
#                  "gene_phase","gene_ID","gene_id")
# colnames(gene_M) <- gene_rename
# full_annotation <- cbind(exon,gene_M)
# col_select<-c(9,1,4,5,6,7,8,10,12,21,16,17,18,19,20,11)
# fa <- full_annotation[,col_select]
# fa$exonSize <-full2$end-full2$start 
# fa$geneSize <-full2$gene_end-full2$gene_start
# exon_place <- NULL
# for(i in 1:length(unique(fa$gene_id))){exon_place <- append(exon_place,c(1:nrow(fa[fa$gene_id == unique(fa$gene_id)[i],])))}
# fa$exon_place <- exon_place
# fa <- fa[,c(1:4,17,19,5:12,18,13:16)]
# rownames(fa) <- fa$exon_ID
# 
# sT_f <- data.frame(sumTab_unstranded,exon_ID=NA,gene_chr=as.character("NA"),
#                    exon_start=as.numeric("NA"),exon_end=as.numeric("NA"),
#                    exonSize=as.numeric("NA"),exon_place=as.numeric("NA"),
#                    score=as.character("NA"),strand=as.character("NA"),
#                    phase=as.character("NA"),Parent=as.character("NA"),
#                    transcript_id=as.character("NA"),gene_ID=as.character("NA"),
#                    gene_start=as.numeric("NA"),gene_end=as.numeric("NA"),
#                    geneSize=as.numeric("NA"),gene_score=as.character("NA"),
#                    gene_strand=as.character("NA"),gene_phase=as.character("NA"),
#                    gene_id=as.character("NA"))
# sT_f[,c(6,7,12:17,21:24)]<-as.character("NA")


#### Creating CpG site subsets based on gene or exon regions ####

# This is currently slightly outdated. I created the `meta_gene` to identify CpGs that intersected with genes. Cataloging
# which gene each genic CpG belonged too (Similar strategy was applied for exons to create `meta_exon`). However, I did this with 
# the stranded data directly from bismark (the cytosine summary flag), but after conversationwith with Steven, it seems that it
# would be better to do this with unstranded CpG data (combining both strand coverage counts). The current approach will like be
# conservative (removing more sites due to insufficient coverage), but potentially reduce (or conversely enhance) the affect of PCR
# biases.

# sprintf("Matching with genes....")
# # Cytosines within genes - variable thresholds
# mC_gene <- mC[meta_gene$cg_index,]
# umC_gene <- umC[meta_gene$cg_index,]
# tC_gene <- tC[meta_gene$cg_index,]
# beta_gene <- beta[meta_gene$cg_index,]
# # Prune by variabel min count thresholds
# tC_gene_min <- apply(tC_gene,1,min) # vector of minimum counts for each cpg
# # At least 1 per sample
# meta_gene_0 <- meta_gene[tC_gene_min>0,]
# mC_gene_0 <- mC_gene[tC_gene_min>0,]
# umC_gene_0 <- umC_gene[tC_gene_min>0,]
# tC_gene_0 <- tC_gene[tC_gene_min>0,]
# beta_gene_0 <- beta_gene[tC_gene_min>0,]
# # At least 5 per sample
# meta_gene_5 <- meta_gene[tC_gene_min>=5,]
# mC_gene_5 <- mC_gene[tC_gene_min>=5,]
# umC_gene_5 <- umC_gene[tC_gene_min>=5,]
# tC_gene_5 <- tC_gene[tC_gene_min>=5,]
# beta_gene_5 <- beta_gene[tC_gene_min>=5,]
# # At least 10 per sample
# meta_gene_10 <- meta_gene[tC_gene_min>=10,]
# mC_gene_10 <- mC_gene[tC_gene_min>=10,]
# umC_gene_10 <- umC_gene[tC_gene_min>=10,]
# tC_gene_10 <- tC_gene[tC_gene_min>=10,]
# beta_gene_10 <- beta_gene[tC_gene_min>=10,]
# 
# saveRDS(mC_gene_5,"Final_mC_gene_5.RData")
# saveRDS(umC_gene_5,"Final_umC_gene_5.RData")
# saveRDS(tC_gene_5,"Final_tC_gene_5.RData")
# saveRDS(beta_gene_5,"Final_beta_gene_5.RData")
# saveRDS(mC_gene_10,"Final_mC_gene_10.RData")
# saveRDS(umC_gene_10,"Final_umC_gene_10.RData")
# saveRDS(tC_gene_10,"Final_tC_gene_10.RData")
# saveRDS(beta_gene_10,"Final_beta_gene_10.RData")
# 
# sprintf("Matching with exons ...")
# # Cytosines with Exon - variable thresholds
# 
# mC_exon <- mC[meta_exon$cg_index,]
# umC_exon <- umC[meta_exon$cg_index,]
# tC_exon <- tC[meta_exon$cg_index,]
# beta_exon <- beta[meta_exon$cg_index,]
# # Prune by variabel min count thresholds
# tC_exon_min <- apply(tC_exon,1,min) # vector of minimum counts for each cpg
# # At least 1 per sample
# meta_exon_0 <- meta_exon[tC_exon_min>0,]
# mC_exon_0 <- mC_exon[tC_exon_min>0,]
# umC_exon_0 <- umC_exon[tC_exon_min>0,]
# tC_exon_0 <- tC_exon[tC_exon_min>0,]
# beta_exon_0 <- beta_exon[tC_exon_min>0,]
# # At least 5 per sample
# meta_exon_5 <- meta_exon[tC_exon_min>=5,]
# mC_exon_5 <- mC_exon[tC_exon_min>=5,]
# umC_exon_5 <- umC_exon[tC_exon_min>=5,]
# tC_exon_5 <- tC_exon[tC_exon_min>=5,]
# beta_exon_5 <- beta_exon[tC_exon_min>=5,]
# # At least 10 per sample
# meta_exon_10 <- meta_exon[tC_exon_min>=10,]
# mC_exon_10 <- mC_exon[tC_exon_min>=10,]
# umC_exon_10 <- umC_exon[tC_exon_min>=10,]
# tC_exon_10 <- tC_exon[tC_exon_min>=10,]
# beta_exon_10 <- beta_exon[tC_exon_min>=10,]
# 
# saveRDS(mC_exon_5,"Final_mC_exon_5.RData")
# saveRDS(umC_exon_5,"Final_umC_exon_5.RData")
# saveRDS(tC_exon_5,"Final_tC_exon_5.RData")
# saveRDS(beta_exon_5,"Final_beta_exon_5.RData")
# saveRDS(mC_exon_10,"Final_mC_exon_10.RData")
# saveRDS(umC_exon_10,"Final_umC_exon_10.RData")
# saveRDS(tC_exon_10,"Final_tC_exon_10.RData")
# saveRDS(beta_exon_10,"Final_beta_exon_10.RData")

## NOTE: This reads in a count matrix from the bedgraph option
 # which by default only considers CpGs

sprintf("Creating summaries...")
#### Summaries ####

# Reading in Data saved from the earlier gene subsetting step
# Currently only reading in info for genes with min coverage of 5 for all samples
setwd("/home/downeyam/Github/2017OAExp_Oysters/input_files/DNAm")
meta_gene_5<-readRDS("Final_meta_gene_5.RData")
mC_gene_5<-readRDS("Final_mC_gene_5.RData")
mC_gene_5 <- data.frame(mC_gene_5)
umC_gene_5<-readRDS("Final_umC_gene_5.RData")
umC_gene_5 <- data.frame(umC_gene_5)
tC_gene_5<-readRDS("Final_tC_gene_5.RData")
tC_gene_5 <- data.frame(tC_gene_5)
beta_gene_5 <- readRDS("Final_beta_gene_5.RData")
beta_gene_5_row <- rownames(beta_gene_5)
beta_gene_5 <- data.frame(beta_gene_5)
rownames(beta_gene_5) <- beta_gene_5_row
# Reading in metadata to subset data later
sample_meta<-readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/meta/metadata_20190811.RData")
sm <- sample_meta[sample_meta$ID != "17099",]
sm_index<-1:nrow(sm)

dim(meta_gene_5)
sum(is.na(meta_gene_5))
# Check to ensure there aren't any NAs
dim(beta_gene_5)
# Same number of loci

# Confirming the rows (CpGs) all match and are sorted the same.
identical(rownames(beta_gene_5),rownames(meta_gene_5))

# Descriptive Methylation 
# # All CpGs
# beta_mean <- rowMeans(beta)
# beta_sd <- unlist(array(beta,1,sd))
# # All CpGs with counts
# beta_0_mean <- rowMeans(beta_0)
# beta_0_sd <- unlist(array(beta_0,1,sd))
# beta_0_CV <- beta_0_sd/beta_0_mean
# beta_5_mean <- rowMeans(beta_5)
# beta_5_sd <- unlist(array(beta_5,1,sd))
# beta_5_CV <- beta_5_sd/beta_5_mean
# beta_10_mean <- rowMeans(beta_10)
# beta_10_sd <- unlist(array(beta_10,1,sd))
# beta_10_CV <- beta_10_sd/beta_10_mean
### CpGs in genes
#beta_gene_0_mean <- rowMeans(beta_gene_0)
#beta_gene_0_sd <- unlist(array(beta_gene_0,1,sd))
#beta_gene_0_CV <- beta_gene_0_sd/beta_gene_0_mean
beta_gene_5_mean <- rowMeans(beta_gene_5)
beta_gene_5_sd <- unlist(array(beta_gene_5,1,sd))
beta_gene_5_CV <- beta_gene_5_sd/beta_gene_5_mean
#beta_gene_10_mean <- rowMeans(beta_gene_10)
#beta_gene_10_sd <- unlist(array(beta_gene_10,1,sd))
#beta_gene_10_CV <- beta_gene_10_sd/beta_gene_10_mean
#### CpG w/ 10threshold - split by condition and time ###


beta_gene_5_control_9 <- beta_gene_5[,sm_index[sm$SFV == "09.400"]]
beta_gene_5_control_80 <- beta_gene_5[,c(sm_index[sm$SFV == "80.400"])]
beta_gene_5_exposed_9 <- beta_gene_5[,sm_index[sm$SFV == "09.2800"]]
beta_gene_5_exposed_80 <- beta_gene_5[,sm_index[sm$SFV == "80.2800"]]
head(beta_gene_5_control_9)
#### CpG w/ 10 threshold sumary split by condition and time ###
# Control - Day 9
beta_gene_5_control_9_mean <- rowMeans(beta_gene_5_control_9)
beta_gene_5_control_9_sd <- unlist(array(beta_gene_5_control_9,1,sd))
beta_gene_5_control_9_CV <- beta_gene_5_control_9_sd/beta_gene_5_control_9_mean
# Control - Day 80
beta_gene_5_control_80_mean <- rowMeans(beta_gene_5_control_80)
beta_gene_5_control_80_sd <- unlist(array(beta_gene_5_control_80,1,sd))
beta_gene_5_control_80_CV <- beta_gene_5_control_80_sd/beta_gene_5_control_80_mean
# Exposed - Day 9
beta_gene_5_exposed_9_mean <- rowMeans(beta_gene_5_exposed_9)
beta_gene_5_exposed_9_sd <- unlist(array(beta_gene_5_exposed_9,1,sd))
beta_gene_5_exposed_9_CV <- beta_gene_5_exposed_9_sd/beta_gene_5_exposed_9_mean
# Exposed - Day 80
beta_gene_5_exposed_80_mean <- rowMeans(beta_gene_5_exposed_80)
beta_gene_5_exposed_80_sd <- unlist(array(beta_gene_5_exposed_80,1,sd))
beta_gene_5_exposed_80_CV <- beta_gene_5_exposed_80_sd/beta_gene_5_exposed_80_mean

# All CpG Gene Summary Table
gene_summary <- data.frame(Gene = meta_gene_5$gene_ID,gene_id=meta_gene_5$gene_id,
                      chr=meta_gene_5$chr_1,pos= meta_gene_5$pos,
                      start=meta_gene_5$start,end=meta_gene_5$end,
                      beta_gene_5_mean,beta_gene_5_sd,beta_gene_5_CV,
                      beta_gene_5_control_9_mean,beta_gene_5_control_9_sd,beta_gene_5_control_9_CV,
                      beta_gene_5_control_80_mean,beta_gene_5_control_80_sd,beta_gene_5_control_80_CV,
                      beta_gene_5_exposed_9_mean,beta_gene_5_exposed_9_sd,beta_gene_5_exposed_9_CV,
                      beta_gene_5_exposed_80_mean,beta_gene_5_exposed_80_sd,beta_gene_5_exposed_80_CV)
#saveRDS(gene_summary,"/home/downeyam/Github/2017OAExp_Oysters/results/DNAm/Final_beta_gene_5_summary.RData")

gl <- data.frame(gene_id=meta_gene_5$gene_id,
                         beta_gene_5_mean,beta_gene_5_sd,beta_gene_5_CV,
                         beta_gene_5_control_9_mean,beta_gene_5_control_9_sd,beta_gene_5_control_9_CV,
                         beta_gene_5_control_80_mean,beta_gene_5_control_80_sd,beta_gene_5_control_80_CV,
                         beta_gene_5_exposed_9_mean,beta_gene_5_exposed_9_sd,beta_gene_5_exposed_9_CV,
                         beta_gene_5_exposed_80_mean,beta_gene_5_exposed_80_sd,beta_gene_5_exposed_80_CV)

test <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/Exon_GeneLoc.RData")
### NEED TO UPDATE FOR EXONS #####

### CpGs in exons
beta_exon_0_mean <- rowMeans(beta_exon_0)
beta_exon_0_sd <- unlist(array(beta_exon_0,1,sd))
beta_exon_0_CV <- beta_exon_0_sd/beta_exon_0_mean
beta_exon_5_mean <- rowMeans(beta_exon_5)
beta_exon_5_sd <- unlist(array(beta_exon_5,1,sd))
beta_exon_5_CV <- beta_exon_5_sd/beta_exon_5_mean
beta_exon_10_mean <- rowMeans(beta_exon_10)
beta_exon_10_sd <- unlist(array(beta_exon_10,1,sd))
beta_exon_10_CV <- beta_exon_10_sd/beta_exon_10_mean
## CpG w/ 10threshold - split by condition and time
beta_exon_10_control_9 <- beta_exon_10[,c(3,5,6,17,18,19)]
beta_exon_10_control_80 <- beta_exon_10[,c(1,2,8,11,22,23)]
beta_exon_10_exposed_9 <- beta_exon_10[,c(7,9,12,13,15,21)]
beta_exon_10_exposed_80 <- beta_exon_10[,c(4,10,14,16,20,24)]
# CpG w/ 10 threshold sumary split by condition and time 
# Control - Day 9
beta_exon_10_control_9_mean <- rowMeans(beta_exon_10_control_9)
beta_exon_10_control_9_sd <- unlist(array(beta_exon_10_control_9,1,sd))
beta_exon_10_control_9_CV <- beta_exon_10_control_9_sd/beta_exon_10_control_9_mean
# Control - Day 80
beta_exon_10_control_80_mean <- rowMeans(beta_exon_10_control_80)
beta_exon_10_control_80_sd <- unlist(array(beta_exon_10_control_80,1,sd))
beta_exon_10_control_80_CV <- beta_exon_10_control_80_sd/beta_exon_10_control_80_mean
# Exposed - Day 9
beta_exon_10_exposed_9_mean <- rowMeans(beta_exon_10_exposed_9)
beta_exon_10_exposed_9_sd <- unlist(array(beta_exon_10_exposed_9,1,sd))
beta_exon_10_exposed_9_CV <- beta_exon_10_exposed_9_sd/beta_exon_10_exposed_9_mean
# Exposed - Day 80
beta_exon_10_exposed_80_mean <- rowMeans(beta_exon_10_exposed_80)
beta_exon_10_exposed_80_sd <- unlist(array(beta_exon_10_exposed_80,1,sd))
beta_exon_10_exposed_80_CV <- beta_exon_10_exposed_80_sd/beta_exon_10_exposed_80_mean

# Gene Summary Table
exon_summary <- cbind(beta_exon_10_mean,beta_exon_10_sd,beta_exon_10_CV,
                      beta_exon_10_control_9_mean,beta_exon_10_control_9_sd,beta_exon_10_control_9_CV,
                      beta_exon_10_control_80_mean,beta_exon_10_control_80_sd,beta_exon_10_control_80_CV,
                      beta_exon_10_exposed_9_mean,beta_exon_10_exposed_9_sd,beta_exon_10_exposed_9_CV,
                      beta_exon_10_exposed_80_mean,beta_exon_10_exposed_80_sd,beta_exon_10_exposed_80_CV)
saveRDS(exon_summary,"Final_beta_exon_summary.RData")

springtf("Creating summary table...")
#### Summary Table ####

nam <- c("Total_CpG",
         "CpG_0","CpG_5x","CpG_10x",
         "CpG_5x_gene","CpG_10x_gene",
         "CpG_5x_exon","CpG_10x_exon")

counts<-c(nrow(tC),
          nrow(tC_0),nrow(tC_5),nrow(tC_10),
          nrow(tC_gene_5),nrow(tC_gene_10),
          nrow(tC_exon_5),nrow(tC_exon_10))

mean_prop_methylation <- c(0.0,
                           mean(beta_0_mean),mean(beta_5_mean),mean(beta_10_mean),
                           mean(beta_gene_5_mean),mean(beta_gene_10_mean),
                           mean(beta_exon_5_mean),mean(beta_exon_10_mean))
sd_prop_methylation <- c(0.0,
                         sd(beta_0_mean),sd(beta_5_mean),sd(beta_10_mean),
                         sd(beta_gene_5_mean),sd(beta_gene_10_mean),
                         sd(beta_exon_5_mean),sd(beta_exon_10_mean))

png("Full_summaryTable.png")
counts_table <- cbind(nam,counts,mean_prop_methylation,sd_prop_methylation)
p <- kableExtra::kable(counts_table) %>% kableExtra::kable_styling()
grid.arrange(p)
dev.off()

sprintf("Printing figures...")
#### Figures ####
## Density plots for different thresholds
# With all CpGs
comb_plot_density <- data.frame(rbind(cbind(beta=beta_mean,Threshold="All"),
                           cbind(beta=beta_0_mean,Threshold="All : Coverage > 0"),
                           cbind(beta=beta_5_mean,Threshold="All : Coverage >= 5"),
                           cbind(beta=beta_10_mean,Threshold="All : Coverage >= 10"),
                           cbind(beta=beta_gene_0_mean,Threshold="Genes : Coverage > 0"),
                           cbind(beta=beta_gene_5_mean,Threshold="Genes : Coverage >= 5"),
                           cbind(beta=beta_gene_10_mean,Threshold="Genes : Coverage >= 10")
))
comb_plot_density$beta <- as.numeric(as.character(comb_plot_density$beta))

ggplot(comb_plot_density,aes(beta,colour=Threshold)) + 
  geom_density() + labs(title = "Mean Methylation Prop")
ggsave("Mean_Methyl_Prop_VariousThresholds.png")


# With just genes
#Means
comb_genes_Mean <- data.frame(rbind(cbind(beta=beta_gene_0_mean,Threshold="Genes : Coverage > 0"),
                                      cbind(beta=beta_gene_5_mean,Threshold="Genes : Coverage >= 5"),
                                      cbind(beta=beta_gene_10_mean,Threshold="Genes : Coverage >= 10")
))
comb_genes_Mean$beta <- as.numeric(as.character(comb_genes_Mean$beta))

ggplot(comb_genes_Mean,aes(beta,colour=Threshold)) + 
  geom_density() + xlim(0,1) +
  labs(title = "Mean Methylation Prop for just genes",xlab="Mean Methylation") 
ggsave("Mean_Methyl_Prop_Genes.png")

#SD
comb_genes_SD <- data.frame(rbind(cbind(beta=beta_gene_0_sd,Threshold="Genes : Coverage > 0"),
                                    cbind(beta=beta_gene_5_sd,Threshold="Genes : Coverage >= 5"),
                                    cbind(beta=beta_gene_10_sd,Threshold="Genes : Coverage >= 10")
))
comb_genes_SD$beta <- as.numeric(as.character(comb_genes_SD$beta))
ggplot(comb_genes_SD,aes(beta,colour=Threshold)) + 
  geom_density() + labs(title = "SD Methylation Prop for just genes") + xlim(0,1)
ggsave("SD_Mean_Methyl_Prop_Genes.png")

## Genes by treatment x time combination
comb_genesByTreatment_Mean <- data.frame(rbind(cbind(beta=beta_gene_10_control_9_mean,Threshold="Control : 9"),
                                  cbind(beta=beta_gene_10_control_80_mean,Threshold="Control : 80"),
                                  cbind(beta=beta_gene_10_exposed_9_mean,Threshold="Exposed : 9"),
                                  cbind(beta=beta_gene_10_exposed_80_mean,Threshold="Exposed : 80"))
)
comb_genesByTreatment_Mean$beta <- as.numeric(as.character(comb_genesByTreatment_Mean$beta))

ggplot(comb_genesByTreatment_Mean,aes(beta,colour=Threshold)) + 
  geom_density() + xlim(0,1) +
  labs(title = "Mean Methylation Prop for just genes by treatment x time (threshold >= 10)",xlab="Mean Methylation") 
ggsave("Mean_Methyl_Prop_Genes_DiffTreatments.png")

## Bar plot of genic methylation levels

comb_GlobalMean <- data.frame(rbind(cbind(beta=mean(beta_gene_10_control_9_mean),sd=sd(beta_gene_10_control_9_mean),Threshold="Control : 09"),
                                    cbind(beta=mean(beta_gene_10_control_80_mean),sd=sd(beta_gene_10_control_80_mean),Threshold="Control : 80"),
                                    cbind(beta=mean(beta_gene_10_exposed_9_mean),sd=sd(beta_gene_10_exposed_9_mean),Threshold="Exposed : 09"),
                                    cbind(beta=mean(beta_gene_10_exposed_80_mean),sd=sd(beta_gene_10_exposed_80_mean),Threshold="Exposed : 80"),
                                    cbind(beta=mean(beta_gene_10_mean),sd=sd(beta_gene_10_mean),Threshold=" All Genes "),
                                    cbind(beta=mean(beta_10_mean),sd=sd(beta_10_mean),Threshold="All CpGs ")
                                               ))
  comb_GlobalMean$beta <- as.numeric(as.character(comb_GlobalMean$beta))
comb_GlobalMean$sd <- as.numeric(as.character(comb_GlobalMean$sd))
ggplot(comb_GlobalMean,aes(x=Threshold,y=beta)) + 
  geom_errorbar(aes(ymin=beta-sd, ymax=beta+sd), width=.1) +
  geom_point() + 
  labs(title = "Global Mean Methylation Prop (threshold >= 10) ") + 
  ylim(0,1)
ggsave("Mean_Methyl_Prop_GlobalMean.png")


