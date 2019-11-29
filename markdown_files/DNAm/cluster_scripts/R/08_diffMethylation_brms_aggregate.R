## Script for aggregating BRMS bayesian binomial model script outputs

# Set the working directory to the folder with the brms scripts outputs you wish to aggregate
setwd("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/processed_samples/07_brmsSummary")
# Create list of each file created by the BRMS script
m_marg <- list.files(pattern="modelMarginalEffects")
m_param <- list.files(pattern="modelParam")
m_sum <- list.files(pattern="modelSummary")
pc <- list.files(pattern="plannedComparisons")

# Subset out range information
# This only need to be done with one of the file lists since the range data is redundant
start <- sub(".*_BRMS_(.+?)_(.+?)_.*","\\1",m_marg,perl=TRUE)
end <- sub(".*_BRMS_(.+?)_(.+?)_.*","\\2",m_marg,perl=TRUE)
# Dataframe of positions
pos <- data.frame(start=as.integer(start),end=as.integer(end))
# Reorder positions
pos <- pos[order(pos$start),]


### Quality Control ###
# Checking to make sure all values exist and do no overlap
all_values <- list()
for(i in 1:nrow(pos)){
  all_values <- list(all_values,seq(pos[i,1],pos[i,2],by = 1))
}
val_list<- unlist(all_values)
if(sum(duplicated(val_list))>0){
  print("Duplicate values detected. List of duplicates:")
  print(val_list[duplicated(val_list)])
  print("This will cause data to be overwritten while combining data outputs")
}

# Create full possible range of values
pos_range <- c(pos[1,1]:pos[nrow(pos),2])
# Check if the ranges in the BRMS output files include all values or if part of the full range
# is missing.
if(sum(is.na(match(pos_range,val_list)))>0){
  print("Range values not continuous. Missing values:")
  print(pos_range[which(is.na(match(pos_range,val_list)))])
}

### Read in first file ###
### The BRMS script was written to generate equal size outputs for each range so 
# be reading in the first dataset, future subsets will be added to these first files.

mme <- readRDS(list.files(pattern=paste0(pos[1,1],"_",pos[1,2],"_modelMarginalEffects")))
ms <- readRDS(list.files(pattern=paste0(pos[1,1],"_",pos[1,2],"_modelSummary")))
mpc <- readRDS(list.files(pattern=paste0(pos[1,1],"_",pos[1,2],"_plannedComparisons")))


for(i in 2:nrow(pos)){
  print(paste0("Processing position range: ",i," of ",nrow(pos)))
  t1<-readRDS(list.files(pattern=paste0(pos[i,1],"_",pos[i,2],"_modelMarginalEffects")))
  t2 <- readRDS(list.files(pattern=paste0(pos[i,1],"_",pos[i,2],"_modelSummary")))
  t3 <- readRDS(list.files(pattern=paste0(pos[i,1],"_",pos[i,2],"_plannedComparisons")))
  
  # Updating marginal effects
  mme$marg_estimate[c(pos[i,1]:pos[i,2]),] <- t1$marg_estimate[c(pos[i,1]:pos[i,2]),]
  mme$se_estimate[c(pos[i,1]:pos[i,2]),] <- t1$se_estimate[c(pos[i,1]:pos[i,2]),]
  mme$l_estimate[c(pos[i,1]:pos[i,2]),] <- t1$l_estimate[c(pos[i,1]:pos[i,2]),]
  mme$u_estimate[c(pos[i,1]:pos[i,2]),] <- t1$u_estimate[c(pos[i,1]:pos[i,2]),]
  
  # Updating model Summary
  ms$Estimate[c(pos[i,1]:pos[i,2]),] <- t2$Estimate[c(pos[i,1]:pos[i,2]),]
  ms$Est.Error[c(pos[i,1]:pos[i,2]),] <- t2$Est.Error[c(pos[i,1]:pos[i,2]),]
  ms$CI.Lower[c(pos[i,1]:pos[i,2]),] <- t2$CI.Lower[c(pos[i,1]:pos[i,2]),]
  ms$CI.Upper[c(pos[i,1]:pos[i,2]),] <- t2$CI.Upper[c(pos[i,1]:pos[i,2]),]
  ms$Rhat[c(pos[i,1]:pos[i,2]),] <- t2$Rhat[c(pos[i,1]:pos[i,2]),]
  ms$Bulk_ESS[c(pos[i,1]:pos[i,2]),] <- t2$Bulk_ESS[c(pos[i,1]:pos[i,2]),]
  ms$Tail_ESS[c(pos[i,1]:pos[i,2]),] <- t2$Tail_ESS[c(pos[i,1]:pos[i,2]),]
  
  # Planned Comparisons
  mpc$Estimate[c(pos[i,1]:pos[i,2]),] <- t3$Estimate[c(pos[i,1]:pos[i,2]),]
  mpc$Est.Error[c(pos[i,1]:pos[i,2]),] <- t3$Est.Error[c(pos[i,1]:pos[i,2]),]
  mpc$CI.Lower[c(pos[i,1]:pos[i,2]),] <- t3$CI.Lower[c(pos[i,1]:pos[i,2]),]
  mpc$CI.Upper[c(pos[i,1]:pos[i,2]),] <- t3$CI.Upper[c(pos[i,1]:pos[i,2]),]  
  mpc$Evid.Ratio[c(pos[i,1]:pos[i,2]),] <- t3$Evid.Ratio[c(pos[i,1]:pos[i,2]),]
  mpc$Post.Prob[c(pos[i,1]:pos[i,2]),] <- t3$Post.Prob[c(pos[i,1]:pos[i,2]),]
  mpc$Significant[c(pos[i,1]:pos[i,2]),] <- t3$Significant[c(pos[i,1]:pos[i,2]),]
}

### Saving Files
saveRDS(mme,"Final_DNAm_gene_BRMS_modelMarginalEffects.RData")
saveRDS(ms,"Final_DNAm_gene_BRMS_modelSummary.RData")
saveRDS(mpc,"Final_DNAm_gene_BRMS_plannedComparisons.RData")



