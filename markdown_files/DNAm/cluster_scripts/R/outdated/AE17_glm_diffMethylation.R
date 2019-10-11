
  library(lme4)
  library(lmerTest)
  library(multcomp)
  library(dplyr)
  library(corrplot)
  
  #install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")
  #devtools::install_github("jakyzhu/MACAU2")
  #install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
  library(MACAU2)
  
  setwd("/home/downeyam/Github/2017OAExp_Oysters/input_files")
  beta <- readRDS("DNAm/Final_beta_gene_10.RData")
  meta <- readRDS("DNAm/Final_meta_gene_10.RData")
  meta_sample <- readRDS("meta/metadata_20190811.RData")
  meta_sample <- meta_sample[meta_sample$ID != "17099",]
  tC <- readRDS("DNAm/Final_tC_gene_5.RData")
  mC <- readRDS("DNAm/Final_mC_gene_10.RData")
  uC <- readRDS("DNAm/Final_umC_gene_5.RData")
  uC <- uC[which(!is.na(uC[,1])),]

meta_sample$tank <- as.factor(meta_sample$tank)
meta_sample$shelf <- as.factor(meta_sample$shelf)
minT_10 <- minT[1:10]
temp <- data.frame(m=unlist(mC[100,]),u=unlist(uC[100,]),
                   SFV=meta_sample$SFV,
                   Pop=meta_sample$Pop,
                   tankID=as.factor(meta_sample$tankID))#,
                   #tank=as.factor(meta_sample$tank),
                   #shelf=as.factor(meta_sample$shelf))

model <- glmer(cbind(m,u)~SFV+(1|tankID),data=temp,family=binomial)
anova(model)
summary(model)
qqnorm(resid(model))
qqline(resid(model))
boxplot(resid(model)~meta_sample$Pop)
boxplot(resid(model)~meta_sample$tankID)
temp$SFV <- as.numeric(temp$SFV)
temp$Pop <- as.numeric(temp$Pop)
temp$tankID <- as.numeric(temp$tankID)
temp$resid <-resid(model)
M <- cor(temp)
corrplot(M)
plot(resid(model)~temp$m)
abline(h=0)
plot(resid(model)~temp$u)
abline(h=0)
temp$m/(temp$m+temp$u)
#### Model Testing and Binomial Regression Exploration with DNAm data ####
  
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
  
#### Testing with MACAU #####

mfit <- macau2(RawCountDataSet=mC[1:10],Phenotypes=meta_sample$treatment,
               LibSize=tC[1:10,],fit.model="BMM",numCore=1)
    
  
#### Analysis ####
diffMeth_pValue <- readRDS("DNAm/DNAm_gene_GLM_modelSummary_wRandTank.RData")
diffMeth_betas <- readRDS("DNAm/Final_beta_gene_5.RData")

plot(diffMeth_pValue$t_score[abs(diffMeth_pValue$t_score[,5]) < 100,5]~diffMeth_pValue$pval[abs(diffMeth_pValue$t_score[,5]) < 100,5])

#Give row names a number for ordering
rownames(diffMeth_pValue$t_score)<-1:nrow(diffMeth_pValue$t_score)
#Reorder tvalues from largest to smallest (absolute values)
diffMeth_pValue$t_score<-diffMeth_pValue$t_score[rev(order(abs(diffMeth_pValue$t_score[,5]))),]
# Create vector of indexs ordered by there descending T values
minT <- as.numeric(rownames(diffMeth_pValue$t_score))
# Loop Through and produce some diagnostic plots looking at the relationship between 
# Methylation (represented as either beta or total counts (tC))
loci <- 1:10
for(j in 1:length(loci)){
  i <- loci[j]
  temp <- data.frame(beta=unlist(diffMeth_betas[minT[i],]),tC=unlist(tC[minT[i],]),
                     SFV=meta_sample$SFV,Pop=meta_sample$Pop,tankID=as.factor(meta_sample$tankID))
  par(mfrow=c(1,1))
  boxplot(tC~SFV,data=temp)
  par(mfrow=c(1,3))
  boxplot(beta~SFV,data=temp,ylim=c(0,1))
  boxplot(beta~Pop,data=temp,ylim=c(0,1))
  boxplot(beta~tankID,data=temp,ylim=c(0,1))
  par(mfrow=c(1,1))
  temp$SFV <- as.numeric(temp$SFV)
  temp$Pop <- as.numeric(temp$Pop)
  temp$tankID <- as.numeric(temp$tankID)
  M <- cor(temp)
  corrplot.mixed(M)
}

## Exploring pvalues after corrections
col<-c("Intercept","SFV80.400","SFV09.2800","SFV80.2800",
       "80.400v09.400","09.2800v09.400","80.2800v09.400",
       "09.2800v80.400","80.2800v80.400","80.2800v09.2800")
colnames(diffMeth_pValue_adj) <-  col
max(diffMeth_pValue$pval[,2])
hist(diffMeth_pValue$pval[,2])
min(diffMeth_pValue_adj$`09.2800v09.400`)
hist(p.adjust(diffMeth_pValue$pval[,2]),method = "fdr")
pVals_corr_09.2800v09.400 <- p.adjust(diffMeth_pValue$pval[,2])
temp <- pVals_corr_09.2800v09.400[pVals_corr_09.2800v09.400 > 0]
temp[temp<0.05]

#### Examining diff Mehtylation vs. diff Expression

diffGene <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/results/DiffExpression_gene_Limma_voom_summaryObject.RData")
head(diffGene$coefficients)
# CvE_D9
plot(-log(diffGene$p.value[,1])~diffGene$coefficients[,1],
     xlab="Log2Fold Change",ylab="-log(p)")
#temp <- !is.na(match(rownames(gene_fit2_contr),topTran$GENEID))
#sub_top_c <- gene_fit2_contr$coefficients[temp,1]
#sub_top_p <- gene_fit2_contr$p.value[temp,1]
points(-log(sub_top_p)~sub_top_c,col="red",pch=16)

