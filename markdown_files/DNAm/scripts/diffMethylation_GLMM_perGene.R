## Libraries
library(lme4)
library(lmerTest)
library(multcomp)
library(svMisc)
source("/home/downeyam/R/basicR_functions.R")

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

setwd("/home/downeyam/Github/2017OAExp_Oysters/input_files")

meta_sample <- readRDS("meta/metadata_20190811.RData")
meta_sample <- meta_sample[meta_sample$ID != "17099",]

# CpGs enriched in genes
beta <- readRDS("DNAm/Final_beta_gene_10.RData")
meta <- readRDS("DNAm/Final_meta_gene_10.RData")
mC <- readRDS("DNAm/Final_mC_gene_10.RData")
mC <- mC[which(!is.na(mC[,1])),]
uC <- readRDS("DNAm/Final_umC_gene_10.RData")
uC <- uC[which(!is.na(uC[,1])),]

genelist <- list()
for( i in 1:length(unique(meta$gene_id))){
  print(paste0("Gene ",i," of ",length(unique(meta$gene_id))))
  Sys.sleep(0.01)
  
  geneID_index <- which(meta$gene_id == unique(meta$gene_id)[i])
  cgIndex <- NULL
  for(j in 1:length(geneID_index)){
    cgIndex <- c(cgIndex,rep(meta$cg_index[geneID_index[j]],nrow(meta_sample)))
  }

  temp_beta<-c(unlist(t(beta[geneID_index,])))
  temp_uC <- c(unlist(t(uC[geneID_index,])))
  temp_mC <- c(unlist(t(mC[geneID_index,])))

  temp <- data.frame(ID=as.character(rep(meta_sample$ID,length(geneID_index))),
                   SFV=rep(meta_sample$SFV,length(geneID_index)),
                   Treatment=rep(meta_sample$Treatment,length(geneID_index)),
                   Time=rep(meta_sample$Time,length(geneID_index)),
                   Pop=rep(meta_sample$Pop,length(geneID_index)),
                   tankID=as.factor(rep(meta_sample$tankID,length(geneID_index))),
                   epf_pH=rep(meta_sample$epf_pH,length(geneID_index)),
                   cg_Index = cgIndex,
                   beta=temp_beta,
                   uC=temp_uC,
                   mC=temp_mC)
  genelist[[i]]<-temp
}

out<-genelist[[2]]
out.glm <- glmer(cbind(mC,uC)~SFV +(1|cg_Index/ID),data=out,family = "binomial")
out.glm <- glmer(beta~SFV +(SFV|cg_Index),data=out,family = "binomial")
out.glm <- glm(beta~SFV,data=out,family = "binomial")
summary(out.glm)
qqnorm(resid(out.glm), main="normal qq-plot, residuals")
qqline(resid(out.glm))
plot(fitted(out.glm), resid(out.glm)) #residuals vs fitted
abline(h=0)

overdisp_fun(out.glm)


# install.packages("R2admb")
# install.packages("glmmADMB", 
#                  repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                          getOption("repos")),
#                  type="source")
# install.packages("betareg")
# library(glmmADMB)
# library(betareg)
# out$cg_Index <- as.factor(out$cg_Index)
# out.beta.glm <- glmmadmb(beta~SFV +(1|cg_Index/ID),data=out,family = "beta",debug = TRUE)

# Scripted loop to apply model to each locus
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
