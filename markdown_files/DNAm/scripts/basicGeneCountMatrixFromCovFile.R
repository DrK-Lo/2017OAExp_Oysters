library(data.table)
library(dplyr)

sprintf("Reading in data and initializing methylation count matrix...")
# Set working directory to base folder where your raw files and bismark outputs reside
setwd("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/singleSample_trimScript_20190719")

# Specific folder where cov matrix are.
bis_outputs <- "bismark/bowtie2"
file_outputs <- "countMatrix"

files <- list.files(bis_outputs)
count_files <- files[grep("deduplicated.bismark.cov.gz",files)]
samples <- substr(count_files,1,5)

init <- fread(paste0(bis_outputs,"/",samples[1],"_DNAm_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"),header=FALSE)
init$comb <- paste0(init$V1,"_",init$V3)

init_C <- cbind(init[,7],init[,5])
colnames(init_C) <- c("comb",paste0(samples[1]))
init_T <- cbind(init[,7],init[,6])
colnames(init_T) <- c("comb",paste0(samples[1]))

sprintf("Starting counting....")
n<- 23
pb <- txtProgressBar(min = 0, max = n, style = 3)

for(i in 2:length(samples)){
  Sys.sleep(0.001)
  setTxtProgressBar(pb, i-1)
  
  temp <- fread(paste0(bis_outputs,"/",samples[i],"_DNAm_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"),header=FALSE)
  temp$comb <- paste0(temp$V1,"_",temp$V3)
  temp_C <- cbind(temp[,7],temp[,5])
  temp_T <- cbind(temp[,7],temp[,6])

  colnames(temp_C) <- c("comb",samples[i])
  colnames(temp_T) <- c("comb",samples[i])
  
  init_C <- left_join(init_C,temp_C,by="comb")
  init_T <- left_join(init_T,temp_T,by="comb")
}
close(pb)

sprintf("Completed merging and reformatting matrices...")
final_C <- data.frame(init_C[,-1])
final_T <- data.frame(init_T[,-1])

row.names(final_C) <- unlist(init_C[,1])
row.names(final_T) <- unlist(init_T[,1])

final_C[is.na(final_C)]<- 0
final_T[is.na(final_T)]<- 0

colnames(final_C) <- samples
colnames(final_T) <- samples

sprintf("Saving files...")
fwrite(final_C,paste0(file_outputs,"/All_methyCountMatrix.csv"))
fwrite(final_T,paste0(file_outputs,"/All_totalCountMatrix.csv"))

saveRDS(final_C,paste0(file_outputs,"/All_methyCountMatrix.RData"))
saveRDS(final_T,paste0(file_outputs,"/All_totalCountMatrix.RData"))


