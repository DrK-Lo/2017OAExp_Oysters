## Libraries
library(lme4)
library(lmerTest)
library(multcomp)
library(svMisc)

## Read in Data
sprintf("Reading in data...")
setwd("/shared_lab/20180226_RNAseq_2017OAExp/DNAm")
beta <- readRDS("trimmedSamples/singleSample_trimScript_20190719/countSummary/Final_beta_gene_5.RData")
meta <- readRDS("trimmedSamples/singleSample_trimScript_20190719/countSummary/Final_meta_gene_5.RData")
meta_sample <- readRDS("metadata/metadata_20190811.RData")
meta_sample <- meta_sample[meta_sample$ID!="17099",]
#tC <- readRDS("trimmedSamples/singleSample_trimScript_20190719/countSummary/Final_tC_gene_5.RData")
#tC <- tC[which(!is.na(tC[,1])),]
mC <- readRDS("trimmedSamples/singleSample_trimScript_20190719/countSummary/Final_mC_gene_5.RData")
mC <- mC[which(!is.na(mC[,1])),]
uC <- readRDS("trimmedSamples/singleSample_trimScript_20190719/countSummary/Final_umC_gene_5.RData")
uC <- uC[which(!is.na(uC[,1])),]


## Performinng logistic regression (glm with binomial family) : Two main explanatory variables  (time and treatment) + interaction
sprintf("Starting model...")
col_p <- c("80.400_09.400","09.2800_09.400","80.2800_09.400","09.2800_80.400","80.2800_80.400","80.2800_09.2800")
#for(i in 1:100){
for(i in 1:nrow(mC)){
  print(paste0("Loci ",i," of ",nrow(mC)))
  Sys.sleep(0.01)
  out <- glm(cbind(unlist(mC[i,]),unlist(uC[i,]))~SFV,data=meta_sample,family = "binomial")
  out.sum <- summary(out)
  ph <- summary(glht(out, linfct = mcp(SFV = "Tukey")))
  if(i == 1){
    model_summary <- list(out.sum)
    model <- list(out)
    est <- cbind(t(unlist(out.sum$coefficients[,1])),t(unlist(ph$test$coefficients)))
    sigma <- cbind(t(unlist(out.sum$coefficients[,2])),t(unlist(ph$test$sigma)))    
    t_score <- cbind(t(unlist(out.sum$coefficients[,3])),t(unlist(ph$test$tstat)))
    p <- cbind(t(unlist(out.sum$coefficients[,4])),t(unlist(ph$test$pvalues)))
   
    }else{
    model_summary <- c(model_summary,list(out.sum))
    model <- c(model,list(out))
    est <- rbind(est,cbind(t(unlist(out.sum$coefficients[,1])),t(unlist(ph$test$coefficients))))
    sigma <- rbind(sigma,cbind(t(unlist(out.sum$coefficients[,2])),t(unlist(ph$test$sigma))))
    t_score <- rbind(t_score,cbind(t(unlist(out.sum$coefficients[,3])),t(unlist(ph$test$tstat))))
    p <- rbind(p,cbind(t(unlist(out.sum$coefficients[,4])),t(unlist(ph$test$pvalues))))
  }
}
model_output <- list(Estimate=est,Sigma=sigma,t_score=t_score,pval=p)

sprintf("Saving outputs...")
saveRDS(model_summary,"DNAm_gene_GLM_modelSummaryObject_fixedEffectsOnly.RData")
saveRDS(model_output,"DNAm_gene_GLM_modelSummary_fixedEffectsOnly.RData")
saveRDS(model,"DNAm_gene_GLM_modelObject_fixedEffectsOnly.RData")
