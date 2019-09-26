
#### Reading in Data a Libraries ####
library(lme4) # mixed effects models
library(multcomp) # post hoc testing
library(devtools)

## Optional : Used for examining zero-inflated negative binomial model
#devtools::install_github("drizopoulos/GLMMadaptive")
#library(GLMMadaptive)
#library(optimx)

setwd("/home/downeyam/Github/2017OAExp_Oysters/input_files")
# This is a list of dataframe with estimates,sigma,t values,
# and p values from running the binomial mixed model 
# (SFV~(1|shelf:tank)) on on all genic CpGs with a coverage of at least 5.
diffMeth_pValue <- readRDS("DNAm/DNAm_gene_GLM_modelSummary_wRandTank.RData")

# Counts Matrices and meta dataframe for all genic CpGs with at least a coverage of 5
mC <- readRDS("DNAm/Final_mC_gene_5.RData")
uC <- readRDS("DNAm/Final_umC_gene_5.RData")
meta <- readRDS("DNAm/Final_meta_gene_5.RData")

# Removing NAs (this should not be necessary, confirm # of rows does not shrink after running)
meta <- meta[which(!is.na(mC[,1])),]
mC <- mC[which(!is.na(mC[,1])),]
uC <- uC[which(!is.na(uC[,1])),]

## Sample data (treatment, site, etc)
meta_sample <- readRDS("meta/metadata_20190811.RData")
meta_sample <- meta_sample[meta_sample$ID!="17099",] # Already removed from count matrices
meta_sample$tank <- as.factor(meta_sample$tank)
meta_sample$shelf <- as.factor(meta_sample$shelf)


#### Calculating Methylation differences between planned comparisons ####
# Creating a per CpG calculation for mean Methylayion difference between planned groups
# Group 1 : Timepoint 9 - 400 vs 2800
# Group 2 : Timepoint 80 - 400 vs 2800
# Group 3 :  All timepoints - 400 vs 2800

# Group 1: Methylation Difference between treatments at day 9
# This will be used for down stream glmer binomial mixed model test
beta<-mC/(mC+uC)
beta.9C <- beta[,meta_sample$SFV == "09.400"]
beta.9E <- beta[,meta_sample$SFV == "09.2800"]
beta.9C_mean <- rowMeans(beta.9C)
beta.9E_mean <- rowMeans(beta.9E)
beta_mean_diff <- beta.9C_mean-beta.9E_mean
beta_mean_diff <- data.frame(index=c(1:length(beta_mean_diff)),beta_diff=beta_mean_diff)
beta_mean_diff_order <- beta_mean_diff[rev(order(beta_mean_diff$beta_diff)),]

# Group 2: Methylation difference between treatments at day 80
beta.80C <- beta[,meta_sample$SFV == "80.400"]
beta.80E <- beta[,meta_sample$SFV == "80.2800"]
beta.80E_mean <- rowMeans(beta.80E)
beta.80C_mean <- rowMeans(beta.80C)
beta_80_mean_diff <- beta.80C_mean-beta.80E_mean
beta_80_mean_diff <- data.frame(index=c(1:length(beta_80_mean_diff)),beta_diff=beta_80_mean_diff)
beta_80_mean_diff_order <- beta_80_mean_diff[rev(order(beta_80_mean_diff$beta_diff)),]

# Group 3: Methylation differences between treatments either day
beta.C <- beta[,meta_sample$treatment == 400]
beta.E <- beta[,meta_sample$treatment == 2800]
beta.E_mean <- rowMeans(beta.E)
beta.C_mean <- rowMeans(beta.C)
beta_all_mean_diff <- beta.C_mean-beta.E_mean
beta_all_mean_diff <- data.frame(index=c(1:length(beta_all_mean_diff)),beta_diff=beta_all_mean_diff)
beta_all_mean_diff_order <- beta_all_mean_diff[rev(order(beta_all_mean_diff$beta_diff)),]
head(beta_all_mean_diff_order)
# Quick visualization of CpG with the largest diff. between treatments
par(mfrow=c(2,2))
locus<-beta_all_mean_diff_order[5,1]
plot(unlist(beta[locus,])~meta_sample$SFV)
boxplot(unlist(beta[locus,])~meta_sample$tankID)
boxplot(unlist(beta[locus,])~meta_sample$Pop)
boxplot(unlist(beta[locus,])~meta_sample$treatment)
meta[locus,]
# Top 3 and 4 are located in same gene  (LOC111124560)
# Top 5 (LOC111106820) - mRNA-cadherin EGF LAG seven-pass G-type receptor 1-like seems interesting

## Using timepoint 9 treatment differences
beta_mean_diff_order[1:10,]
#         index  beta_diff
#4863813   68199 0.6901395
#21405600 250806 0.6542031
#4645597   66443 0.6512425
#16656845 224675 0.6330357
#15451380 217314 0.6161995
#2953601   40923 0.6136095
#19162747 238883 0.6106158
#3988824   59960 0.6033316
#18539674 233871 0.6011208
#4778328   67234 0.6004323

# Keep top differences for glmm
top_09_diff <- c(beta_mean_diff_order$index[1:10])

# Subset count matrices to the top candidates
uC_top <- uC[top_09_diff,] 
mC_top <- mC[top_09_diff,]
# Visualize some of the top candidates
beta_top<-mC_top/(uC_top+mC_top)
par(mfrow=c(2,2))
candidate<-10 # Change this form 1-10 to switch between top candidates
plot(unlist(beta_top[candidate,])~meta_sample$SFV)
boxplot(unlist(beta_top[candidate,])~meta_sample$tankID)
boxplot(unlist(beta_top[candidate,])~meta_sample$Pop)
boxplot(unlist(beta_top[candidate,])~meta_sample$treatment)

# Confirm converage for top candidate is just has no methylation in our 09.2800 treatment
(uC_top+mC_top)[1,] # Looks like it isn't a coverage issue
beta_top[1,meta_sample$SFV == "09.2800"]
beta_top[1,meta_sample$SFV == "09.400"]


#### GLMER Model with permutations ####
# THIS WILL TAKE TIME DEPENDING ON # of PERMUTATIONS
# Saved RData object below was based on 1000 permutations

## Running glmer first on actual data then performing a permutation (1000) to look at the random range of statistics
 # Fstat - F statistic of an anova on the model
 # t - t statistics from the multcomp function based on the planned comparison
 # p - p value from the multcomp function based on the planned comparison

PERM<-1000
sprintf("Starting model...")
col_p <- c("80.400_09.400","09.2800_09.400","80.2800_09.400","09.2800_80.400","80.2800_80.400","80.2800_09.2800")
for(i in 1:length(top_09_diff)){
  tryCatch({
    print(paste0("Loci ",i," of ",length(top_09_diff)))
    Sys.sleep(0.01)
    # Run model using actual data and store all outputs
    out <- glmer(cbind(unlist(mC_top[i,]),unlist(uC_top[i,]))~SFV+(1|tank:shelf),data=meta_sample,family = "binomial")
    #out <- glmer.nb(cbind(unlist(mC_top[1,]),unlist(uC_top[1,]))~SFV+(1|tank:shelf),data=meta_sample)
    out.sum <- summary(out)
    ph <- summary(glht(out, linfct = mcp(SFV = "Tukey")))
    if(i == 1){
      model_summary <- list(out.sum)
      model <- list(out)
      est <- cbind(t(unlist(out.sum$coefficients[,1])),t(unlist(ph$test$coefficients)))
      sigma <- cbind(t(unlist(out.sum$coefficients[,2])),t(unlist(ph$test$sigma)))
      t_score <- cbind(t(unlist(out.sum$coefficients[,3])),t(unlist(ph$test$tstat)))
      p <- cbind(t(unlist(out.sum$coefficients[,4])),t(unlist(ph$test$pvalues)))
      Fstat <- anova(out)[1,4] 
    }else{
      model_summary <- c(model_summary,list(out.sum))
      model <- c(model,list(out))
      est <- rbind(est,cbind(t(unlist(out.sum$coefficients[,1])),t(unlist(ph$test$coefficients))))
      sigma <- rbind(sigma,cbind(t(unlist(out.sum$coefficients[,2])),t(unlist(ph$test$sigma))))
      t_score <- rbind(t_score,cbind(t(unlist(out.sum$coefficients[,3])),t(unlist(ph$test$tstat))))
      p <- rbind(p,cbind(t(unlist(out.sum$coefficients[,4])),t(unlist(ph$test$pvalues))))
      Fstat <- c(Fstat,anova(out)[1,4])
    }
    # Perform permutations by reshuffling response variables (methyl and unmethylated counts)
    for(j in 1:PERM){
      print(paste0("Perm ",j," of ",PERM))
      Sys.sleep(0.01)
      mC_temp <- sample(mC_top[1,],size = ncol(mC_top),replace = FALSE)
      uC_temp <- sample(uC_top[1,],size = ncol(uC_top),replace = FALSE)
      out <- glmer(cbind(unlist(mC_temp),unlist(uC_temp))~SFV+(1|tank:shelf),data=meta_sample,family = "binomial")
      #out <- glmer.nb(cbind(unlist(mC_top[1,]),unlist(uC_top[1,]))~SFV+(1|tank:shelf),data=meta_sample)
      out.sum <- summary(out)
      ph <- summary(glht(out, linfct = mcp(SFV = "Tukey")))
      
      if(i == 1){
        if(j == 1){
          Fstat_perm <- matrix(0,ncol=PERM,nrow=length(top_09_diff))
          t_score_perm <- matrix(0,ncol=PERM,nrow=length(top_09_diff))
          p_perm <- matrix(0,ncol=PERM,nrow=length(top_09_diff))
        }
        Fstat_perm[1,j] <- anova(out)[1,4]
        t_score_perm[1,j] <- unlist(ph$test$tstat[2])
        p_perm[1,j] <- unlist(ph$test$pvalues[2])
      } else{
        Fstat_perm[i,j] <- anova(out)[1,4]
        t_score_perm[i,j] <- unlist(ph$test$tstat[2])
        p_perm[i,j] <- unlist(ph$test$pvalues[2])
        }
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

### Save outputs
# perm_GLMM_outputs <- list(model_summary,model,est,sigma,t_score,p,Fstat,Fstat_perm,t_score_perm,p_perm)
# saveRDS(perm_GLMM_outputs,"DNAm/perm1000_test_GLMM_topDiff_09vTreatments.RData")

### Read in Data
perm_GLMM_outputs <- readRDS("DNAm/perm1000_test_GLMM_topDiff_09vTreatments.RData")  

## Note on object
 # Its a list with both lists and matrices
 # [[1]] $model_summary : A list of the model summaries (i.e. summary())  for actual data for each locus
 # [[2]] $model : A list of the model objects generated by running the glmer function with actual data for each locus
 # [[3]] $est : Matrix of estimates for the 3 primary estimates and follow up planned comparisons for actual data for each locus
 # [[4]] $sigma : Matrix of estimates for the 3 primary sigmas and follow up planned comparisons for actual data for each locus
 # [[5]] $tvalue : Matrix of estimates for the 3 primary t scores and follow up planned comparisons for actual data for each locus
 # [[6]] $p : Matrix of estimates for the 3 primary p values and follow up planned comparisons for actual data for each locus
 # [[7]] $Fstat : Vector of F statistics from an anova using the actual data for each locus
 # [[8]] $Fstat_perm : Matrix of f stats based on reshuffling actual response data (rows= # of permuted values, col= number of cpgs)
 # [[9]] $t_score_perm : Matrix of t scores based on reshuffling actual response data (rows= # of permuted values, col= number of cpgs)
 # [[10]] $p_perm : Matrix of p values based on reshuffling actual response data (rows= # of permuted values, col= number of cpgs)

#### Summarizing outputs #### 

### Code below will need to be tweaked if you are starting from the saved perm_GLMM_output object
 ## It assumes you are running strait from the model.
 ## You will need to specific the object names + the specific list item (e.g. perm_GLMM_outputs$Fstat_perm)

# Look at histograms of the F, t, and p values distribution
par(mfrow=c(1,3))
perm_hist<-1
hist(Fstat_perm[perm_hist,])
hist(t_score_perm[perm_hist,])
hist(p_perm[perm_hist,])

# Calculating the quantile for each statistic based on permutation
custom_quant <- function(x) {quantile(x,type=1,probs = c(0,0.01,0.05,0.25,0.5,0.75,0.95,0.99,1))}
F_range <- t(apply(Fstat_perm,1,function(x){custom_quant(x)}))
t_range <- t(apply(t_score_perm,1,function(x){custom_quant(x)}))
p_range <- t(apply(p_perm,1,function(x){custom_quant(x)}))
# Finalizing permutation summary tables by adding the actual statistic value on the right

# Relabel and add tentative significance column
col <-  c("q_0","q_0.01","q_0.05","q_0.25","q_0.5","q_0.75","q_0.95","q_0.99","q_1","Actual")
Fsummary <- data.frame(F_range,Actual_F=Fstat)
colnames(Fsummary) <-  col
Fsummary$Significant <- Fsummary$Actual>Fsummary$q_0.95
Tsummary <- data.frame(t_range,Actual_t=t_score[,6])
colnames(Tsummary) <-  col
Tsummary$Significant <- (Tsummary$Actual>Tsummary$q_0.95 | Tsummary$Actual<Tsummary$q_0.05)
psummary <- data.frame(p_range,Actual_p=p[,6])
colnames(psummary) <-  col

# Print out summarized tables
Fsummary
Tsummary
psummary

##### Single locus examination with multiple models ####

# create temp new data frame with both response (counts) and explanatory (trtTime factor and tankID) variables
# Unaltered
temp_data <- data.frame(uc= unlist(uC_top[1,]),mc = unlist(mC_top[1,]),beta=unlist(beta_top[1,]),
                        TrtTime = meta_sample$SFV,tankID=meta_sample$tankID)
# Adding methylated C counts to add variation into the '09.2800' treatment methylation (current all inds are 0% methylated)
temp_data_2 <- data.frame(uc= unlist(uC_top[1,]),mc = unlist(mC_top[1,]),beta=unlist(beta_top[1,]),
                          TrtTime = meta_sample$SFV,tankID=meta_sample$tankID)
temp_data_2$mc[temp_data_2$TrtTime=="09.2800"] <- c(0,0,0,0,0,1)
## Thoughts : We tested adding minimal variation into our % methylatation estimate within a TimexTreatmet and 
 # this seems to solve the issue with convergence (no error appeared when running glmer with rand effects).
 # This also seemed important in the context of the the general binomial model (with or without rand effects)
 # which wasn't able to appropriate handle a factor level with no variation within the response variable. 
 # Even adding a small amount of variation (adding a single methylated cytosine to one individual), was sufficient
 # to improve this issue.


### Different models
# Binomial Regression w/ random effects (default) 
glmer_count <- glmer(cbind(uc,mc) ~ TrtTime + (1 | tankID), data = temp_data_2,
                     family = binomial)
# Issues with non-convergence

# Basic Binomial with no random effects 
glm_count <- glm(cbind(uc,mc) ~ TrtTime, data = temp_data_2,
                     family = binomial)
summary(glm_count)
ph <- summary(glht(glm_count, linfct = mcp(TrtTime = "Tukey")))
ph$test$tstat
ph$test$pvalues

# Binomial Regression w/ random effects (optimx controller) supposedly helps with
# non-convergence
glmer_count <- glmer(cbind(uc,mc) ~ TrtTime + (1 | tankID), data = temp_data,
                     family = binomial,
                     glmerControl(optimizer ='optimx', 
                                  optCtrl=list(method='nlminb')))

# Convert binomial counts into single response (beta) and running negative binomial
glmer_beta <- glmer(beta ~ TrtTime + (1 | tankID), data = temp_data,
                  family = negative.binomial)

# Using package glmmAdaptive to explore negative binomial models and zero inflation
nb <- mixed_model(fixed = beta ~ TrtTime, random = ~ 1 | tankID, data = temp_data,
                  family = negative.binomial())
nb_zi <- mixed_model(fixed = beta ~ TrtTime, random = ~ 1 | tankID, data = temp_data,
                  family = zi.negative.binomial(),zi_fixed = ~1)
summary(nb_zi)


