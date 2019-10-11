library(data.table)
library(dplyr)
library(matrixStats)

sprintf("Reading in data and initializing methylation count matrix...")
# Set working directory to base folder where your raw files and bismark outputs reside
setwd("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/processed_samples")

# Specific folder where cov matrix are.
input_folder  <- "03_CytoSummaries"
output_folder <- "04_countMatrix"

files <- list.files(input_folder)
count_files <- files[grep("CytoSummary.CpG_report.txt",files)]
samples <- substr(count_files,1,5)

init <- fread(paste0(input_folder,"/",samples[1],"_CytoSummary.CpG_report.txt"),header=FALSE)

# This columns correspond to the CpG metadata. This should only need to be stored once
# since it will be the same for all samples
# 1 - Chromosome
# 2 - position
# 3 - strand
# 6 - Cytosine motif
# 7 - Trinucleotide motif
sum.table <- init[,c(1,2,3,6,7)]
row_nam <- paste0(sum.table$V1,"_",sum.table$V2)

# Going with the initialization step
init_C <- matrix(0,ncol=length(samples),nrow=nrow(init))
colnames(init_C) <- samples
rownames(init_C) <- row_nam
init_U <- matrix(0,ncol=length(samples),nrow=nrow(init))
colnames(init_U) <- samples
rownames(init_U) <- row_nam

sprintf("Starting counting....")
n<- length(samples)
pb <- txtProgressBar(min = 0, max = n, style = 3)

for(i in 1:length(samples)){
  Sys.sleep(0.001)
  setTxtProgressBar(pb, i)
  
  temp <- fread(paste0(input_folder,"/",samples[i],"_CytoSummary.CpG_report.txt"),header=FALSE)
  
  init_C[,i] <- temp$V4
  init_U[,i] <- temp$V5
}
close(pb)

sprintf("Completed merging and reformatting matrices...")

# Calculating the total count coverage for each cytosine
init_T <- init_C + init_U
colnames(init_T) <- samples
rownames(init_T) <- row_nam
# Calculating individual level prop. methylation (# methylated/total cytosines)
init_B<- init_C/init_T
init_B[is.na(init_B)]<-0
colnames(init_B) <- samples
rownames(init_B) <- row_nam

## Destranding data and summing coverage between paried cytosines within a CpG.

# Separate by positive and negative strand
sum.table.pos <- sum.table[sum.table$V3 == "+",]
sum.table.neg <- sum.table[sum.table$V3 == "-",]
# Create a consolidated destranded metadata object
sum.table.destrand <- data.frame(chr=sum.table.pos$V1,
                                 start=sum.table.pos$V2,end=sum.table.neg$V2,
                                 C_motif=sum.table.pos$V6,Tri_motif=sum.table.pos$V7)
# Check the cytosines between strands are match up
diff <- abs(sum.table.destrand$end-sum.table.destrand$start)
sum(diff>1)
# If this is 0 then they are matching properly

# Destrand and sum cytosine between strands
init_U_destrand <- init_U[sum.table$V3 == "+",] + init_U[sum.table$V3 == "-",]
init_C_destrand <- init_C[sum.table$V3 == "+",] + init_C[sum.table$V3 == "-",]
init_T_destrand <- init_T[sum.table$V3 == "+",] + init_T[sum.table$V3 == "-",]
init_B_destrand <- init_C_destrand/init_T_destrand
init_B_destrand[is.na(init_B_destrand)]<-0

# ## Histogram just to evaluate the coverage equality between strands
# # Lets only look at loci that meet the 5 min for all samples coverage critera
# T_destrand_rowMin <- rowMins(init_T_destrand)
# # Count number that meet the minimum critera (5 coverage for each sample)
# sum(T_destrand_rowMin>=5) 
# # Determine prop of counts for a CpG that are from positive strand
# div <- (init_T[sum.table$V3 == "+",]/init_T_destrand)*100
# div[is.na(div)]<-0
# div <- div[T_destrand_rowMin >= 5,]
# div_mean<- rowMeans(div)
# div_sd <- rowSds(div)
# T_count <- rowSums(init_T_destrand[T_destrand_rowMin >= 5,])
# B_prop <- rowMeans(init_B_destrand[T_destrand_rowMin >= 5,])
# # Create 2x2 panel plot visualizing any potential bias between strands with summary stats
# png(paste0(output_folder,"positiveStrandContribution_Hist.png"))
# par(mfrow=c(2,2))
# hist(div,xlim=c(0.01,100),xlab="% from Postive Strand",
#      main="% positive strand contribution to destranded count")
# hist(div_mean,xlim=c(0.01,100),xlab="Mean % from Postive Strand (for all samples)",
#      main="Mean % positive strand across samples")
# hist(div_sd,xlim=c(0.01,100),xlab="SD of % from Postive Strand (for all samples)",
#      main="SD of % positive strand across samples")
# plot(div_sd~T_count,xlab="Summed CpG Coverage",ylab="SD of % positive strand across samples")
# dev.off()
# # Single panel plot visualizing potential bias on the estimated methylation level (B)
# png(paste0(output_folder,"positiveStrandEffectOnBeta.png"))
# png("positiveStrandEffectOnBeta.png")
# par(mfrow=c(1,1))
# plot(B_prop~div_mean,xlab="Mean % positive strand across samples",main="Mean CpG beta vs Mean % positive strand")
# dev.off()

sprintf("Saving files...")
## Stranded
# Summary Table
fwrite(sum.table,paste0(output_folder,"/CG_stranded_summaryTable.csv"))
saveRDS(sum.table,paste0(output_folder,"/CG_stranded_summaryTable.RData"))
# CSVs
fwrite(init_C,paste0(output_folder,"/CG_stranded_methylCountMatrix.csv"))
fwrite(init_U,paste0(output_folder,"/CG_stranded_unmethylCountMatrix.csv"))
fwrite(init_T,paste0(output_folder,"/CG_stranded_TotalCountMatrix.csv"))
fwrite(init_B,paste0(output_folder,"/CG_stranded_BetaMatrix.csv"))
# RData
saveRDS(init_C,paste0(output_folder,"/CG_stranded_methylCountMatrix.RData"))
saveRDS(init_U,paste0(output_folder,"/CG_stranded_unmethylCountMatrix.RData"))
saveRDS(init_T,paste0(output_folder,"/CG_stranded_TotalCountMatrix.RData"))
saveRDS(init_B,paste0(output_folder,"/CG_stranded_BetaMatrix.RData"))

## Unstranded
# Summary Table
fwrite(sum.table.destrand,paste0(output_folder,"/CG_unstranded_summaryTable.csv"))
saveRDS(sum.table.destrand,paste0(output_folder,"/CG_unstranded_summaryTable.RData"))
# CSVs
fwrite(init_C_destrand,paste0(output_folder,"/CG_unstranded_methylCountMatrix.csv"))
fwrite(init_U_destrand,paste0(output_folder,"/CG_unstranded_unmethylCountMatrix.csv"))
fwrite(init_T_destrand,paste0(output_folder,"/CG_unstranded_TotalCountMatrix.csv"))
fwrite(init_B_destrand,paste0(output_folder,"/CG_unstranded_BetaMatrix.csv"))
# RData
saveRDS(init_C_destrand,paste0(output_folder,"/CG_unstranded_methylCountMatrix.RData"))
saveRDS(init_U_destrand,paste0(output_folder,"/CG_unstranded_unmethylCountMatrix.RData"))
saveRDS(init_T_destrand,paste0(output_folder,"/CG_unstranded_TotalCountMatrix.RData"))
saveRDS(init_B_destrand,paste0(output_folder,"/CG_unstranded_BetaMatrix.RData"))

