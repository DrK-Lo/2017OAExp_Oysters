
  library(lme4)
  library(lmerTest)
  library(multcomp)
  library(dplyr)
  setwd("/home/downeyam/Github/2017OAExp_Oysters/input_files")
  beta <- readRDS("DNAm/Final_beta_gene_10.RData")
  meta <- readRDS("DNAm/Final_meta_gene_10.RData")
  meta_sample <- readRDS("meta/metadata_20190811.RData")
  meta_sample <- meta_sample[meta_sample$ID != "17099",]
  tC <- readRDS("DNAm/Final_tC_gene_10.RData")
  mC <- readRDS("DNAm/Final_mC_gene_10.RData")
  uC <- readRDS("DNAm/Final_umC_gene_5.RData")
  uC <- uC[which(!is.na(uC[,1])),]

beta2 <- data.frame(ID=colnames(beta),beta=t(beta[1,]))
sample_test <- full_join(beta2,meta_sample)
# sample_test <- data.frame(meta_sample$sample_name,beta=t(beta[1,]),
#                             meta_sample$Time,meta_sample$Treatment,meta_sample$shelf,meta_sample$tank,meta_sample$tankID,
#                             meta_sample$Pop)
boxplot(sample_test$beta~sample_test$Treatment*sample_test$shelf)
boxplot(sample_test$beta~sample_test$Treatment*sample_test$tank)
boxplot(sample_test$beta~sample_test$Treatment*sample_test$tankID,las=2) 
boxplot(sample_test$beta~sample_test$Treatment*sample_test$Pop*sample_test$tankID,las=2) 
boxplot(sample_test$beta~sample_test$Time*sample_test$Pop,las=2) 
# Some differences at least in first CpG with random effects
table(meta_sample$tankID,meta_sample$Pop)
table(meta_sample$Treatment,meta_sample$Pop)
table(meta_sample$Time,meta_sample$Pop)
table(meta_sample$tankID,meta_sample$timepoint)
# Mostly Balanced
Full <- glmer(cbind(unlist(mC[1,]),unlist(uC[1,]))~Treatment*Time+(Treatment|Pop),data=meta_sample,family = "binomial")
summary(Full)
# Initial Exploration of different GLM models
Full <- glmer(cbind(unlist(mC[1,]),unlist(uC[1,]))~Treatment*Time+(1|shelf/tank)+(1|Pop),data=meta_sample,family = "binomial")
summary(Full)
qqnorm(Full)
## These next two models are equivalent (1|tankID) vs. (1|tank:shelf)
Full <- glmer(cbind(unlist(mC[1,]),unlist(uC[1,]))~Treatment*Time+(1|tankID)+(1|Pop),data=meta_sample,family = "binomial")
summary(Full)
Full <- glmer(cbind(unlist(mC[1,]),unlist(uC[1,]))~Treatment*Time+(1|tank:shelf)+(1|Pop),data=meta_sample,family = "binomial")
summary(Full)
qqnorm(resid(Full), main="normal qq-plot, residuals")
qqline(resid(Full))
plot(fitted(Full), resid(Full)) #residuals vs fitted
abline(h=0)

# Just pop as random effect
Full2 <- glmer(cbind(unlist(mC[1,]),unlist(uC[1,]))~Treatment*Time+(1|Pop),data=meta_sample,family = "binomial")
summary(Full2)
qqnorm(resid(Full2), main="normal qq-plot, residuals")
qqline(resid(Full2))
plot(fitted(Full2), resid(Full)) #residuals vs fitted
abline(h=0)
# Just the fixed effects
Partial <- glm(cbind(unlist(mC[1,]),unlist(tC[1,]- mC[1,]))~Treatment*Time,data=meta_sample,family = "binomial")
Partial <- glm(cbind(unlist(mC[1,]),unlist(tC[1,]- mC[1,]))~SFV,data=meta_sample,family = "binomial")
summary(Partial)
qqnorm(resid(Partial), main="normal qq-plot, residuals")
qqline(resid(Partial))
plot(fitted(Partial), resid(Full)) #residuals vs fitted
abline(h=0)
ph <- summary(glht(Partial, linfct = mcp(SFV = "Tukey")))
ph
# Examine residuals of simple (partail model) based on chapter 5 (page 150) of Zuur book
  
# Summaries
  out.sum <- summary(Partial)
  ph <- summary(glht(Partial, linfct = mcp(SFV = "Tukey")))
# Estimate
  est <- cbind(t(unlist(out.sum$coefficients[,1])),t(unlist(ph$test$coefficients)))
  est <- rbind(est,cbind(t(unlist(out.sum$coefficients[,1])),t(unlist(ph$test$coefficients))))
# Sigma
  sigma <- cbind(t(unlist(out.sum$coefficients[,2])),t(unlist(ph$test$sigma)))
  sigma <- rbind(sigma,cbind(t(unlist(out.sum$coefficients[,2])),t(unlist(ph$test$sigma))))
# Tstat (or z score)
  t_score <- cbind(t(unlist(out.sum$coefficients[,3])),t(unlist(ph$test$tstat)))
  t_score <- rbind(t_score,cbind(t(unlist(out.sum$coefficients[,3])),t(unlist(ph$test$tstat))))
# Tstat (or z score)
  p <- cbind(t(unlist(out.sum$coefficients[,4])),t(unlist(ph$test$pvalues)))
  p <- rbind(p,cbind(t(unlist(out.sum$coefficients[,4])),t(unlist(ph$test$pvalues))))
# Removing the randoms
  summary(Partial)
  anova(Full,Partial)
  anova(Full,Full2)
  anova(Random_Tank,Partial)

### Current version of script run on server for only fixed effects
  ## Read in Data
  sprintf("Reading in data...")
  setwd("/shared_lab/20180226_RNAseq_2017OAExp/DNAm")
  beta <- readRDS("trimmedSamples/singleSample_trimScript_20190719/countSummary/Final_beta_gene_10.RData")
  meta <- readRDS("trimmedSamples/singleSample_trimScript_20190719/countSummary/Final_meta_gene_10.RData")
  meta_sample <- readRDS("metadata/metadata_20190811.RData")
  tC <- readRDS("trimmedSamples/singleSample_trimScript_20190719/countSummary/Final_tC_gene_10.RData")
  tC <- tC[which(!is.na(tC[,1])),]
  mC <- readRDS("trimmedSamples/singleSample_trimScript_20190719/countSummary/Final_mC_gene_10.RData")
  mC <- mC[which(!is.na(mC[,1])),]
  uC <- readRDS("trimmedSamples/singleSample_trimScript_20190719/countSummary/Final_umC_gene_10.RData")
  uC <- uC[which(!is.na(uC[,1])),]
  
  
  ## Performinng logistic regression (glm with binomial family) : Two main explanatory variables  (time and treatment) + interaction
  sprintf("Starting model...")
  col_p <- c("80.400_09.400","09.2800_09.400","80.2800_09.400","09.2800_80.400","80.2800_80.400","80.2800_09.2800")
  for(i in 1:100){
    #for(i in 1:nrow(tC)){
    print(paste0("Loci ",i," of ",nrow(tC)))
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
  model_output <- list(est,sigma,t_score,p)
  
  sprintf("Saving outputs...")
  saveRDS(model_summary,"DNAm_gene_GLM_modelSummaryObject_fixedEffectsOnly.RData")
  saveRDS(model_output,"DNAm_gene_GLM_modelSummary_fixedEffectsOnly.RData")
  saveRDS(model,"DNAm_gene_GLM_modelObject_fixedEffectsOnly.RData")
  
diffMeth_pValue <- readRDS("DNAm/PlannedComparison_pValueRaw.RData")
diffMeth_pValue <- readRDS("DNAm/DNAm_gene_GLM_modelSummary_fixedEffectsOnly.RData")
head(diffMeth_pValue$pval)
diffMeth_pValue_adj <- data.frame(matrix(0,ncol=ncol(diffMeth_pValue),nrow=nrow(diffMeth_pValue)))
dim(diffMeth_pValue_adj)
diffMeth_pValue$`09.2800_09.400`
## Perform FDR p value adjustment
for(i in 1:ncol(diffMeth_pValue)){
 diffMeth_pValue_adj[,i] <- p.adjust(diffMeth_pValue[,i],method = "fdr")  
}
colnames(diffMeth_pValue_adj) <- colnames(diffMeth_pValue)
hist(diffMeth_pValue_adj$`09.2800_09.400`)  
hist(diffMeth_pValue_adj$`80.2800_80.400`)
hist(diffMeth_pValue_adj$`80.400_09.400`)

#### Examining diff Mehtylation vs. diff Expression

diffGene <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/results/DiffExpression_gene_Limma_voom_summaryObject.RData")

# CvE_D9
plot(-log(diffGene$p.value[,1])~diffGene$coefficients[,1],
     xlab="Log2Fold Change",ylab="-log(p)")
#temp <- !is.na(match(rownames(gene_fit2_contr),topTran$GENEID))
#sub_top_c <- gene_fit2_contr$coefficients[temp,1]
#sub_top_p <- gene_fit2_contr$p.value[temp,1]
3points(-log(sub_top_p)~sub_top_c,col="red",pch=16)

