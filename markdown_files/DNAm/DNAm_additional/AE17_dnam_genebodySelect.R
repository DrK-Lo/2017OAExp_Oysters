library(dplyr)

# Script was used to generate feature specific subsets of the DNA methylation data.


#### Data #### 
# Set working directory
setwd("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/processed_samples")
# Read in the total count matrixes by type (beta, total, summary, etc)
cpg_sum <- readRDS("05_countSummary/CG_unstranded_summaryTable_5.RData")
beta <- readRDS("05_countSummary/CG_unstranded_beta_5.RData")
total <- readRDS("05_countSummary/CG_unstranded_tC_5.RData")
mC <- readRDS("05_countSummary/CG_unstranded_mC_5.RData")
umC <- readRDS("05_countSummary/CG_unstranded_umC_5.RData")
# Annotation file with all cpgs and all correponding features it belongs too.
annot <- read.table("08_geneLevelSummary/DNAm_Data/CpGsAnnotated.bed",
                       header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")


#### Processing Data #####  

## Create unique label for subsetting
# Create unique cpg ids to match with
annot$ID <- paste0(annot$V1,"_",annot$V2)
cpg_sum$ID <- paste0(cpg_sum$chr,"_",cpg_sum$start)

## Left join to combine summary table for each unique cpg and total annotation file.
# Combine annotations 
comb <- left_join(cpg_sum,annot)

## Subsetting 
# Keep only cpgs in genes
comb_gene <- comb[comb$V11 == "gene",]
# Keep only cpgs in exons
comb_exon <- comb[comb$V11 == "exon",]

### NOTE: these will multiple copies (rows) of the same CpG IF
# that cpg was found in multiple genes or (in the case of exons) the exons 
# within a gene overlapped in the annotation file. The later instance seems
# to be how we have more rows in the complete exon file than in the gene file.
comb_gene_unique <- comb_gene[!duplicated(comb_gene$ID),]
comb_exon_unique <- comb_exon[!duplicated(comb_exon$ID),]


## Filter count matrices using the non duplicate ID from the combined data 
## For genes
out <- rownames(beta) %in% comb_gene_unique$ID
beta_gene <- beta[out,]
# tC filter
out <- rownames(total) %in% comb_gene_unique$ID
total_gene <- total[out,]
# mC filter
out <- rownames(mC) %in% comb_gene_unique$ID
mC_gene <- mC[out,]
# umC filter
out <- rownames(umC) %in% comb_gene_unique$ID
umC_gene <- umC[out,]
## For exons
out <- rownames(beta) %in% comb_exon_unique$ID
beta_exon <- beta[out,]
# tC filter
out <- rownames(total) %in% comb_exon_unique$ID
total_exon <- total[out,]
# mC filter
out <- rownames(mC) %in% comb_exon_unique$ID
mC_exon <- mC[out,]
# umC filter
out <- rownames(umC) %in% comb_exon_unique$ID
umC_exon <- umC[out,]

#### Visualizing the data ####

## Gene level
# All beta (for each individual and cpg) included separately on hist
png("05_countSummary/histogram_unstranded_5min_geneBody_All_BetaValues.png")
hist(beta_gene)
dev.off()
# Mean beta for each cpg (averaged across all individuals)
png("05_countSummary/histogram_unstranded_5min_geneBody_Mean_BetaValues.png")
hist(rowMeans(beta_gene))
dev.off()
# Mean beta for each cpg (averaged within treatment)
meta<-readRDS("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/metadata/metadata_20190811.RData")
meta<-meta[meta$ID!="17099",]
C9 <- beta_gene[,meta$SFV=="09.400"]
C80 <- beta_gene[,meta$SFV=="80.400"]
E9 <- beta_gene[,meta$SFV=="09.2800"]
E80 <- beta_gene[,meta$SFV=="80.2800"]
#mean by treatment
C9_mean <- rowMeans(C9)
C80_mean <- rowMeans(C80)
E9_mean <- rowMeans(E9)
E80_mean <- rowMeans(E80)

blue.trans <- rgb(0, 0, 1, 1/3)
red.trans <- rgb(1, 0, 0, 1/3)
#plot
png("05_countSummary/histogram_unstranded_5min_geneBody_byTreatment_Mean_BetaValues.png")
hist(rowMeans(C9),col = blue.trans,breaks=100)
par(new = TRUE)
hist(rowMeans(E9),col = red.trans, axes = F, xlab = "", ylab = "", main = "",breaks=100)
dev.off()
#plot density 
png("05_countSummary/density_unstranded_5min_geneBody_byTreatment_Mean_BetaValues.png",
    res = 400,width = 20,height = 20,units = "cm")
plot(density(rowMeans(C80)),col = col_perm[2],
     main="DNA methylation density by treatment-time levels",
     xlab="Prop. DNA Methylation")
lines(density(rowMeans(E9)),col = col_perm[3])
lines(density(rowMeans(C9)),col = col_perm[1])
lines(density(rowMeans(E80)),col = col_perm[4])
abline(v=median(rowMeans(C9)),col = col_perm[1])
abline(v=median(rowMeans(C80)),col = col_perm[2])
abline(v=median(rowMeans(E9)),col = col_perm[3])
abline(v=median(rowMeans(E80)),col = col_perm[4])
dev.off()

#Stacked barplot 
#color palatte
library("RColorBrewer")
pal <- brewer.pal(n = 12, name = 'Paired')
col_perm <- c(pal[1:2],pal[5:6],pal[12])
#quantile groups
qq_C9<-cut(C9_mean,
        breaks=seq(min(C9_mean),max(C9_mean), length.out=21),
        include.lowest=T, labels=F)
qq_C80<-cut(C80_mean,
           breaks=seq(min(C80_mean),max(C80_mean), length.out=21),
           include.lowest=T, labels=F)
qq_E9<-cut(E9_mean,
           breaks=seq(min(E9_mean),max(E9_mean), length.out=21),
           include.lowest=T, labels=F)
qq_E80<-cut(E80_mean,
           breaks=seq(min(E80_mean),max(E80_mean), length.out=21),
           include.lowest=T, labels=F)
# matrices for storage
means_quantile <- matrix(ncol=4,nrow=length(unique(qq)))

upperCI_quantile <- matrix(ncol=4,nrow=length(unique(qq)))
lowerCI_quantile <- matrix(ncol=4,nrow=length(unique(qq)))
ci <- function(x){1.96*(sd(x)/sqrt(length(x)))}
# Calculate mean methylation by quantile
for(i in 1:length(unique(qq))){
  C9_sum <- sum(C9_mean[qq_C9==i])
  C80_sum <- sum(C80_mean[qq_C80==i])
  E9_sum <- sum(E9_mean[qq_E9==i])
  E80_sum <- sum(E80_mean[qq_E80==i])
  means_quantile[i,] <- c(C9_sum,C80_sum,E9_sum,E80_sum) 
}
#plot
png("05_countSummary/stackedbar_unstranded_5min_geneBody_byTreatment_Mean_BetaValues.png")
barplot(t(means_quantile),beside = TRUE,
        col = c(col_perm[1],col_perm[2],col_perm[3],col_perm[4]),
        xlab="Methylation (percentiles")
dev.off()

## Saving files
#genes
saveRDS(comb_gene,"05_countSummary/CG_unstranded_summaryTable_geneOnly_5.RData")
saveRDS(beta_gene,"05_countSummary/CG_unstranded_beta_geneOnly_5.RData")
saveRDS(total_gene,"05_countSummary/CG_unstranded_tC_geneOnly_5.RData")
saveRDS(mC_gene,"05_countSummary/CG_unstranded_mC_geneOnly_5.RData")
saveRDS(umC_gene,"05_countSummary/CG_unstranded_umC_geneOnly_5.RData")
#exons
saveRDS(comb_exon,"05_countSummary/CG_unstranded_summaryTable_exonOnly_5.RData")
saveRDS(beta_exon,"05_countSummary/CG_unstranded_beta_exonOnly_5.RData")
saveRDS(total_exon,"05_countSummary/CG_unstranded_tC_exonOnly_5.RData")
saveRDS(mC_exon,"05_countSummary/CG_unstranded_mC_exonOnly_5.RData")
saveRDS(umC_exon,"05_countSummary/CG_unstranded_umC_exonOnly_5.RData")

