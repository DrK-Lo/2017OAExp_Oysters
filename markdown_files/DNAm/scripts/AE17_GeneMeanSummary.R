
### Libaries
library(edgeR)
library(dplyr)
library(ggplot2)
library(multcomp)
library(matrixStats)
source("/home/downeyam/R/basicR_functions.R")
### Data
setwd("/home/downeyam/Github/2017OAExp_Oysters/input_files")

meta_sample <- readRDS("meta/metadata_20190811.RData")
meta_sample <- meta_sample[meta_sample$ID != "17099",]
# All CpG enriched (using sum 5 threshold)
beta_total_sum5 <- readRDS("DNAm/Final_beta_sum5.RData")
# All CpGs enriched 
beta_total_0 <- readRDS("DNAm/Final_beta_0.RData")
beta_total_5 <- readRDS("DNAm/Final_beta_5.RData")
beta_total_10 <- readRDS("DNAm/Final_beta_10.RData")
# CpGs enriched in genes
beta_gene_0 <- readRDS("DNAm/Final_beta_gene_0.RData")
beta_gene_5 <- readRDS("DNAm/Final_beta_gene_5.RData")
beta_gene_10 <- readRDS("DNAm/Final_beta_gene_10.RData")
# CpGs enriched in exons
beta_exon_0 <- readRDS("DNAm/Final_beta_exon_0.RData")
beta_exon_0 <- beta_exon_0[which(!is.na(beta_exon_0[,1])),]
beta_exon_5 <- readRDS("DNAm/Final_beta_exon_5.RData")
beta_exon_5 <- beta_exon_5[which(!is.na(beta_exon_5[,1])),]
beta_exon_10 <- readRDS("DNAm/Final_beta_exon_10.RData")
beta_exon_10 <- beta_exon_10[which(!is.na(beta_exon_10[,1])),]
# Meta Data files
meta_10 <- readRDS("DNAm/Final_meta_10.RData")
meta_gene_10 <- readRDS("DNAm/Final_meta_gene_10.RData")
meta_exon_10 <- readRDS("DNAm/Final_meta_exon_10.RData")
meta_exon_10 <- meta_exon_10[which(!is.na(meta_exon_10$exon_index)),]

beta_perLoci_summary <- rbind(
data.frame(beta=beta_total_0,label="Total_0"),
data.frame(beta=beta_total_0,label="Total_0"),
data.frame(beta=beta_total_0,label="Total_0")
)
# Checks to make sure dataframe matches up
setdiff(rownames(meta_exon_10),rownames(beta_exon_10))

# Subseting introns
beta_intron_0 <- beta_gene_0[setdiff(rownames(beta_gene_0),rownames(beta_exon_0)),]
beta_intron_5 <- beta_gene_5[setdiff(rownames(beta_gene_5),rownames(beta_exon_5)),]
beta_intron_10 <- beta_gene_10[setdiff(rownames(beta_gene_10),rownames(beta_exon_10)),]

# RNA matrix after min cpm = 1 per gene per sample threshold filtering
rna <- readRDS("RNA/STAR_mapping/RSEM_output/GeneCountMatrix_Filtered_DGEListObj.RData")
rna_c <- rna$counts
rna_cpm <- cpm(rna_c)

# Considering intergenic (gene) methylation for loci with >10 coverage per sample
## Global Methylation levels by treate x time
# All CpGs
C9_beta <- subset(beta_total_10,select=which(meta_sample$SFV == "09.400"))
C80_beta <- subset(beta_total_10,select=which(meta_sample$SFV == "80.400"))
E9_beta <- subset(beta_total_10,select=which(meta_sample$SFV == "09.2800"))
E80_beta <- subset(beta_total_10,select=which(meta_sample$SFV == "80.2800"))
# CpGs in Genes
C9_gene_beta <- subset(beta_gene_10,select=which(meta_sample$SFV == "09.400"))
C80_gene_beta <- subset(beta_gene_10,select=which(meta_sample$SFV == "80.400"))
E9_gene_beta <- subset(beta_gene_10,select=which(meta_sample$SFV == "09.2800"))
E80_gene_beta <- subset(beta_gene_10,select=which(meta_sample$SFV == "80.2800"))
# CpGs in Exons
C9_exon_beta <- subset(beta_exon_10,select=which(meta_sample$SFV == "09.400"))
C80_exon_beta <- subset(beta_exon_10,select=which(meta_sample$SFV == "80.400"))
E9_exon_beta <- subset(beta_exon_10,select=which(meta_sample$SFV == "09.2800"))
E80_exon_beta <- subset(beta_exon_10,select=which(meta_sample$SFV == "80.2800"))
# CpGs in Introns
C9_intron_beta <- subset(beta_intron_10,select=which(meta_sample$SFV == "09.400"))
C80_intron_beta <- subset(beta_intron_10,select=which(meta_sample$SFV == "80.400"))
E9_intron_beta <- subset(beta_intron_10,select=which(meta_sample$SFV == "09.2800"))
E80_intron_beta <- subset(beta_intron_10,select=which(meta_sample$SFV == "80.2800"))

### Means by for each sample 
# Dist. Total Methylation hist.
#Total_beta_10 <- rowMeans(as.matrix(beta_total_10),na.rm=TRUE)
#hist(Total_beta_10)
# Total methylation
C9_beta_mean <- colMeans(as.matrix(C9_beta),na.rm = TRUE)
C80_beta_mean <- colMeans(as.matrix(C80_beta),na.rm = TRUE)
E9_beta_mean <- colMeans(as.matrix(E9_beta),na.rm = TRUE)
E80_beta_mean <- colMeans(as.matrix(E80_beta),na.rm = TRUE)
# Gene methylation
C9_gene_beta_mean <- colMeans(as.matrix(C9_gene_beta),na.rm = TRUE)
C80_gene_beta_mean <- colMeans(as.matrix(C80_gene_beta),na.rm = TRUE)
E9_gene_beta_mean <- colMeans(as.matrix(E9_gene_beta),na.rm = TRUE)
E80_gene_beta_mean <- colMeans(as.matrix(E80_gene_beta),na.rm = TRUE)
# Exon methylation
C9_exon_beta_mean <- colMeans(as.matrix(C9_exon_beta),na.rm = TRUE)
C80_exon_beta_mean <- colMeans(as.matrix(C80_exon_beta),na.rm = TRUE)
E9_exon_beta_mean <- colMeans(as.matrix(E9_exon_beta),na.rm = TRUE)
E80_exon_beta_mean <- colMeans(as.matrix(E80_exon_beta),na.rm = TRUE)
# Intron methylation
C9_intron_beta_mean <- colMeans(as.matrix(C9_intron_beta),na.rm = TRUE)
C80_intron_beta_mean <- colMeans(as.matrix(C80_intron_beta),na.rm = TRUE)
E9_intron_beta_mean <- colMeans(as.matrix(E9_intron_beta),na.rm = TRUE)
E80_intron_beta_mean <- colMeans(as.matrix(E80_intron_beta),na.rm = TRUE)

##Summary Stats
# Total Genome
beta_sampleGlobalMean <- data.frame(beta=c(C9_beta_mean,C80_beta_mean,E9_beta_mean,E80_beta_mean),
           Treatment=c(rep("Ambient",11),rep("OA 2800",12)),
           Time=c(rep("09",6),rep("80",5),rep("09",6),rep("80",6)),
           SFV = as.factor(c(rep("C9",6),rep("C80",5),rep("E9",6),rep("E80",6))),
           Feature="All")
g_beta_lm <- lm(beta~Treatment*Time,data=beta_sampleGlobalMean)
anova(g_beta_lm)
g_beta_lm_s <- lm(beta~SFV,data=beta_sampleGlobalMean)
anova(g_beta_lm_s)
summary(glht(g_beta_lm_s, linfct = mcp(SFV = "Tukey")))
ggplot(beta_sampleGlobalMean,aes(y=beta,group=interaction(Treatment,Time),colour=Treatment)) + geom_boxplot()
# Gene
beta_gene_sampleGlobalMean <- data.frame(beta=c(C9_gene_beta_mean,C80_gene_beta_mean,E9_gene_beta_mean,E80_gene_beta_mean),
                                    Treatment=c(rep("Ambient",11),rep("OA 2800",12)),
                                    Time=c(rep("09",6),rep("80",5),rep("09",6),rep("80",6)),
                                    SFV = as.factor(c(rep("C9",6),rep("C80",5),rep("E9",6),rep("E80",6))),
                                    Feature="Gene")
g_gene_beta_lm <- lm(beta~Treatment*Time,data=beta_gene_sampleGlobalMean)
anova(g_gene_beta_lm)
g_gene_beta_lm_s <- lm(beta~SFV,data=beta_gene_sampleGlobalMean)
anova(g_gene_beta_lm_s)
summary(glht(g_gene_beta_lm_s, linfct = mcp(SFV = "Tukey")))
ggplot(beta_gene_sampleGlobalMean,aes(y=beta,group=interaction(Treatment,Time),colour=Treatment)) + geom_boxplot()
# Exon
beta_exon_sampleGlobalMean <- data.frame(beta=c(C9_exon_beta_mean,C80_exon_beta_mean,E9_exon_beta_mean,E80_exon_beta_mean),
                                         Treatment=c(rep("Ambient",11),rep("OA 2800",12)),
                                         Time=c(rep("09",6),rep("80",5),rep("09",6),rep("80",6)),
                                         SFV = as.factor(c(rep("C9",6),rep("C80",5),rep("E9",6),rep("E80",6))),
                                         Feature="Exon")
g_exon_beta_lm <- lm(beta~Treatment*Time,data=beta_exon_sampleGlobalMean)
anova(g_exon_beta_lm)
g_exon_beta_lm_s <- lm(beta~SFV,data=beta_exon_sampleGlobalMean)
anova(g_exon_beta_lm_s)
summary(glht(g_exon_beta_lm_s, linfct = mcp(SFV = "Tukey")))
ggplot(beta_exon_sampleGlobalMean,aes(y=beta,group=interaction(Treatment,Time),colour=Treatment)) + geom_boxplot()
# Intron
beta_intron_sampleGlobalMean <- data.frame(beta=c(C9_intron_beta_mean,C80_intron_beta_mean,E9_intron_beta_mean,E80_intron_beta_mean),
                                         Treatment=c(rep("Ambient",11),rep("OA 2800",12)),
                                         Time=c(rep("09",6),rep("80",5),rep("09",6),rep("80",6)),
                                         SFV = as.factor(c(rep("C9",6),rep("C80",5),rep("E9",6),rep("E80",6))),
                                         Feature="Intron")
g_intron_beta_lm <- lm(beta~Treatment*Time,data=beta_intron_sampleGlobalMean)
anova(g_intron_beta_lm)
g_intron_beta_lm_s <- lm(beta~SFV,data=beta_intron_sampleGlobalMean)
anova(g_intron_beta_lm_s)
summary(glht(g_intron_beta_lm_s, linfct = mcp(SFV = "Tukey")))
ggplot(beta_intron_sampleGlobalMean,aes(y=beta,group=interaction(Treatment,Time),colour=Treatment)) + geom_boxplot()

beta_comb_sampleGlobalMean <- data.frame(rbind(beta_sampleGlobalMean,
                                               beta_gene_sampleGlobalMean,
                                               beta_exon_sampleGlobalMean,
                                               beta_intron_sampleGlobalMean))
beta_comb_sampleGlobalMean$comb <- as.factor(paste0(beta_comb_sampleGlobalMean$SFV,"_",beta_comb_sampleGlobalMean$Feature))
g_comb_beta_lm <- lm(beta~Treatment*Time*Feature,data=beta_comb_sampleGlobalMean)
anova(g_comb_beta_lm)
g_comb_beta_lm <- lm(beta~comb,data=beta_comb_sampleGlobalMean)
anova(g_comb_beta_lm)
#summary(glht(g_comb_beta_lm, linfct = mcp(comb = "Tukey")))

ggplot(beta_comb_sampleGlobalMean,
       aes(y=beta*100,x=SFV,
           group=interaction(Time,Treatment))) + 
   geom_boxplot(outlier.alpha = 0) + 
  geom_jitter() + 
  facet_grid(.~Feature) +
  theme_bw() +
  scale_x_discrete(limits=c("C9","E9","C80","E80")) +
  labs(x="Treatment.Time",y="Mean Global Methylation (%)",colour="Time.Treatment") +
  theme(strip.text.x = element_text(size=12),
        strip.background = element_rect(colour=NULL, fill="white"))

# Grand mean of intragenic methylation
C9_beta_gmean <- mean(C9_beta_mean)
C80_beta_gmean <- mean(C80_beta_mean)
E9_beta_gmean <- mean(E9_beta_mean)
E80_beta_gmean <- mean(E80_beta_mean)
# se
# Global CI 95%
C9_beta_se <- se(C9_beta_mean)
C80_beta_se <- se(C80_beta_mean) 
E9_beta_se <- se(E9_beta_mean) 
E80_beta_se <- se(E80_beta_mean) 
# Global CI 95%
C9_beta_ci <- se(C9_beta_mean)  * 1.96
C80_beta_ci <- se(C80_beta_mean) * 1.96
E9_beta_ci <- se(E9_beta_mean) * 1.96
E80_beta_ci <- se(E80_beta_mean) * 1.96
# Global Methylation summary table and plot
g_sum <- data.frame(label=c("0400.09","0400.80","2800.09","2800.80"),
                    mean=c(C9_beta_gmean,C80_beta_gmean,E9_beta_gmean,E80_beta_gmean),
                    se=c(C9_beta_se,C80_beta_se,E9_beta_se,E80_beta_se),
                    ci=c(C9_beta_ci,C80_beta_ci,E9_beta_ci,E80_beta_ci))
#saveRDS(g_sum,"/home/downeyam/Github/2017OAExp_Oysters/input_files/DNAm/globalMethylation_gene_10_summary.RData")
ggplot(g_sum,aes(y=mean*100,x=label)) + 
  geom_point(size=5) + 
  geom_errorbar(aes(ymin=(mean-ci)*100,ymax=(mean+ci)*100),width=.2) + 
  scale_x_discrete(limits=c("0400.09","2800.09","0400.80","2800.80")) +
  ylim(78,83) + labs(x="Treatment and Time",y="Mean Global Genic Methlyation (%)",size=6) + theme_bw()

C_beta_gmean <- mean(c(C9_beta_mean,C80_beta_mean))
E_beta_gmean <-  mean(c(E9_beta_mean,E80_beta_mean))
C_beta_ci <- se(c(C9_beta_mean,C80_beta_mean))#*1.96
E_beta_ci <- se(c(E9_beta_mean,E80_beta_mean))#*1.96

g_sum <- data.frame(label=c("Control","Exposed"),
                    mean=c(C_beta_gmean,E_beta_gmean),
                    ci=c(C_beta_ci,E_beta_ci))
ggplot(g_sum,aes(y=mean*100,x=label)) + 
  geom_point(size=5) + 
  geom_errorbar(aes(ymin=(mean-ci)*100,ymax=(mean+ci)*100),width=.2) + 
  #scale_x_discrete(size=10)+
  ylim(78,83) + labs(x="Treatment and Time",y="Mean Global Methlyation (%)",size=6) + theme_bw()


## Coefficient of Variation 
# Sample CV of global methylation
# Total genome
C9_beta_CV <-  apply(as.matrix(C9_beta),1,function(x)sd(x)/mean(x)) 
C80_beta_CV <- apply(as.matrix(C80_beta),1,function(x)sd(x)/mean(x))
E9_beta_CV <- apply(as.matrix(E9_beta),1,function(x)sd(x)/mean(x))
E80_beta_CV <- apply(as.matrix(E80_beta),1,function(x)sd(x)/mean(x))
# Gene
C9_gene_beta_CV <-  apply(as.matrix(C9_gene_beta),1,function(x)sd(x)/mean(x)) 
C80_gene_beta_CV <- apply(as.matrix(C80_gene_beta),1,function(x)sd(x)/mean(x))
E9_gene_beta_CV <- apply(as.matrix(E9_gene_beta),1,function(x)sd(x)/mean(x))
E80_gene_beta_CV <- apply(as.matrix(E80_gene_beta),1,function(x)sd(x)/mean(x))
# Exon
C9_exon_beta_CV <-  apply(as.matrix(C9_exon_beta),1,function(x)sd(x)/mean(x)) 
C80_exon_beta_CV <- apply(as.matrix(C80_exon_beta),1,function(x)sd(x)/mean(x))
E9_exon_beta_CV <- apply(as.matrix(E9_exon_beta),1,function(x)sd(x)/mean(x))
E80_exon_beta_CV <- apply(as.matrix(E80_exon_beta),1,function(x)sd(x)/mean(x))
# Intron
C9_intron_beta_CV <-  apply(as.matrix(C9_intron_beta),1,function(x)sd(x)/mean(x)) 
C80_intron_beta_CV <- apply(as.matrix(C80_intron_beta),1,function(x)sd(x)/mean(x))
E9_intron_beta_CV <- apply(as.matrix(E9_intron_beta),1,function(x)sd(x)/mean(x))
E80_intron_beta_CV <- apply(as.matrix(E80_intron_beta),1,function(x)sd(x)/mean(x))

#Total
C9_dat <- data.frame(beta=C9_beta_CV,Treatment="Ambient",Time="09",SFV="C9",Feature="All")
C80_dat <- data.frame(beta=C80_beta_CV,Treatment="Ambient",Time="80",SFV="C80",Feature="All")
E9_dat <- data.frame(beta=E9_beta_CV,Treatment="OA 2800",Time="09",SFV="E9",Feature="All")
E80_dat <- data.frame(beta=E80_beta_CV,Treatment="OA 2800",Time="80",SFV="E80",Feature="All")

C80_dat$beta <- C80_dat$beta-C9_dat$beta
E9_dat$beta <- E9_dat$beta-C9_dat$beta
E80_dat$beta <- E80_dat$beta-C9_dat$beta
C9_dat$beta <- C9_dat$beta-C9_dat$beta

beta_sampleGlobalCV <- data.frame(rbind(C9_dat,C80_dat,E9_dat,E80_dat))
g_beta_lm <- lm(beta~Treatment*Time,data=beta_sampleGlobalCV)
anova(g_beta_lm)
ggplot(beta_sampleGlobalCV,aes(y=beta,group=interaction(Treatment,Time),colour=Treatment)) + 
  geom_boxplot() + ylim(-0.2,0.2)

# Gene
beta_gene_sampleGlobalCV <- data.frame(beta=c(C9_gene_beta_CV,C80_gene_beta_CV,E9_gene_beta_CV,E80_gene_beta_CV),
                                       Treatment=c(rep("Ambient",11),rep("OA 2800",12)),
                                       Time=c(rep("09",6),rep("80",5),rep("09",6),rep("80",6)),
                                       SFV = as.factor(c(rep("C9",6),rep("C80",5),rep("E9",6),rep("E80",6))),
                                       Feature="Gene")
g_gene_beta_lm <- lm(beta~Treatment*Time,data=beta_gene_sampleGlobalCV)
anova(g_gene_beta_lm)
g_gene_beta_lm_s <- lm(beta~SFV,data=beta_gene_sampleGlobalCV)
anova(g_gene_beta_lm_s)
summary(glht(g_gene_beta_lm_s, linfct = mcp(SFV = "Tukey")))
ggplot(beta_gene_sampleGlobalCV,aes(y=beta,group=interaction(Treatment,Time),colour=Treatment)) + geom_boxplot()
# Original used when CV summary was for each sample NOT each loci
# beta_sampleGlobalCV <- data.frame(beta=c(C9_beta_CV,C80_beta_CV,E9_beta_CV,E80_beta_CV),
#                                   Treatment=c(rep("Ambient",11),rep("OA 2800",12)),
#                                   Time=c(rep("09",6),rep("80",5),rep("09",6),rep("80",6)),
#                                   SFV = as.factor(c(rep("C9",6),rep("C80",5),rep("E9",6),rep("E80",6))),
#                                   Feature="All")
# g_beta_lm <- lm(beta~Treatment*Time,data=beta_sampleGlobalCV)
# anova(g_beta_lm)
# ggplot(beta_sampleGlobalCV,aes(y=beta,group=interaction(Treatment,Time),colour=Treatment)) + geom_boxplot()
# # Gene
# beta_gene_sampleGlobalCV <- data.frame(beta=c(C9_gene_beta_CV,C80_gene_beta_CV,E9_gene_beta_CV,E80_gene_beta_CV),
#                                          Treatment=c(rep("Ambient",11),rep("OA 2800",12)),
#                                          Time=c(rep("09",6),rep("80",5),rep("09",6),rep("80",6)),
#                                          SFV = as.factor(c(rep("C9",6),rep("C80",5),rep("E9",6),rep("E80",6))),
#                                          Feature="Gene")
# g_gene_beta_lm <- lm(beta~Treatment*Time,data=beta_gene_sampleGlobalCV)
# anova(g_gene_beta_lm)
# g_gene_beta_lm_s <- lm(beta~SFV,data=beta_gene_sampleGlobalCV)
# anova(g_gene_beta_lm_s)
# summary(glht(g_gene_beta_lm_s, linfct = mcp(SFV = "Tukey")))
# ggplot(beta_gene_sampleGlobalCV,aes(y=beta,group=interaction(Treatment,Time),colour=Treatment)) + geom_boxplot()
# Exon
beta_exon_sampleGlobalCV <- data.frame(beta=c(C9_exon_beta_CV,C80_exon_beta_CV,E9_exon_beta_CV,E80_exon_beta_CV),
                                         Treatment=c(rep("Ambient",11),rep("OA 2800",12)),
                                         Time=c(rep("09",6),rep("80",5),rep("09",6),rep("80",6)),
                                         SFV = as.factor(c(rep("C9",6),rep("C80",5),rep("E9",6),rep("E80",6))),
                                         Feature="Exon")
g_exon_beta_lm <- lm(beta~Treatment*Time,data=beta_exon_sampleGlobalCV)
anova(g_exon_beta_lm)
g_exon_beta_lm_s <- lm(beta~SFV,data=beta_exon_sampleGlobalCV)
anova(g_exon_beta_lm_s)
summary(glht(g_exon_beta_lm_s, linfct = mcp(SFV = "Tukey")))
ggplot(beta_exon_sampleGlobalCV,aes(y=beta,group=interaction(Treatment,Time),colour=Treatment)) + geom_boxplot()
# Intron
beta_intron_sampleGlobalCV <- data.frame(beta=c(C9_intron_beta_CV,C80_intron_beta_CV,E9_intron_beta_CV,E80_intron_beta_CV),
                                           Treatment=c(rep("Ambient",11),rep("OA 2800",12)),
                                           Time=c(rep("09",6),rep("80",5),rep("09",6),rep("80",6)),
                                           SFV = as.factor(c(rep("C9",6),rep("C80",5),rep("E9",6),rep("E80",6))),
                                           Feature="Intron")
g_intron_beta_lm <- lm(beta~Treatment*Time,data=beta_intron_sampleGlobalCV)
anova(g_intron_beta_lm)
g_intron_beta_lm_s <- lm(beta~SFV,data=beta_intron_sampleGlobalCV)
anova(g_intron_beta_lm_s)
summary(glht(g_intron_beta_lm_s, linfct = mcp(SFV = "Tukey")))
ggplot(beta_intron_sampleGlobalCV,aes(y=beta,group=interaction(Treatment,Time),colour=Treatment)) + geom_boxplot()

beta_comb_sampleGlobalCV <- data.frame(rbind(beta_sampleGlobalCV,
                                               beta_gene_sampleGlobalCV,
                                               beta_exon_sampleGlobalCV,
                                               beta_intron_sampleGlobalCV))
beta_comb_sampleGlobalCV$comb <- as.factor(paste0(beta_comb_sampleGlobalCV$SFV,"_",beta_comb_sampleGlobalCV$Feature))
g_comb_beta_lm <- lm(beta~Treatment*Time*Feature,data=beta_comb_sampleGlobalCV)
anova(g_comb_beta_lm)
g_comb_beta_lm <- lm(beta~comb,data=beta_comb_sampleGlobalCV)
anova(g_comb_beta_lm)
#summary(glht(g_comb_beta_lm, linfct = mcp(comb = "Tukey")))

ggplot(beta_comb_sampleGlobalCV,
       aes(y=beta*100,x=SFV,
           group=interaction(Time,Treatment))) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter() + 
  facet_grid(.~Feature) +
  theme_bw() +
  scale_x_discrete(limits=c("C9","E9","C80","E80")) +
  labs(x="Treatment.Time",y="CV Global Methylation (%)",colour="Time.Treatment") +
  theme(strip.text.x = element_text(size=12),
        strip.background = element_rect(colour=NULL, fill="white"))

# Multiplotting mean and CV 
#Mean

a <- ggplot(beta_comb_sampleGlobalMean,
       aes(y=as.numeric(beta*100),x=SFV,
           group=interaction(Time,Treatment))) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter() + 
  facet_grid(.~Feature) +
  theme_bw() +
  scale_x_discrete(limits=c("C9","E9","C80","E80"),breaks=NULL) +
  labs(x=NULL,y="Mean Global Methylation (%)",colour="Time.Treatment") +
  theme(strip.text.x = element_text(size=12),
        panel.grid.major.x = element_blank() ,
        strip.background = element_rect(colour=NULL, fill="white"))
#CV
b <- ggplot(beta_comb_sampleGlobalCV,
       aes(y=beta,x=SFV,
           group=interaction(Time,Treatment))) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter() + 
  facet_grid(.~Feature) +
  theme_bw() +
  scale_x_discrete(limits=c("C9","E9","C80","E80")) +
  labs(x="Treatment.Time",y="Global Methylation CV",colour="Time.Treatment") +
  theme(strip.text.x = element_text(size=12),
        panel.grid.major.x = element_blank() ,
        strip.background = element_rect(colour=NULL, fill="white"))
library(cowplot)
plot_grid(a, b, ncol=1, align="v")


## Summarizing methylation by gene
rm(DNAm_gene_summary)
#for(i in 1:10){
for(i in 1:length(unique(meta$gene_id))){
  print(paste0("Processing gene...",i," of ",length(unique(meta$gene_id))))
  meta_temp <- subset(meta,meta$gene_id == unique(meta$gene_id)[i])
  beta_temp <- as.matrix(subset(beta,meta$gene_id == unique(meta$gene_id)[i]))

  if(nrow(beta_temp) > 0){
    beta_temp_C09 <- as.matrix(beta_temp[,meta_sample_no99$SFV == "09.400"])
    beta_temp_C80 <- as.matrix(beta_temp[,meta_sample_no99$SFV == "09.2800"])
    beta_temp_E09 <- as.matrix(beta_temp[,meta_sample_no99$SFV == "80.400"])
    beta_temp_E80 <- as.matrix(beta_temp[,meta_sample_no99$SFV == "80.2800"])
    size <- meta_temp$end[1] - meta_temp$start[1]
    CpG_Range <- max(meta_temp$pos)-min(meta_temp$pos)
    Span <- CpG_Range/size
    #All samples
    b_median <- mean(rowMedians(beta_temp))
    b_mean <- mean(rowMeans(beta_temp))
    b_mean_sd <- mean(rowSds(beta_temp))
    b_sd_mean <- sd(rowMeans(beta_temp))
    b_mean_CV <-  mean(rowSds(beta_temp)/rowMeans(beta_temp))
    b_CV_mean <- sd(rowMeans(beta_temp))/mean(rowMeans(beta_temp))
    b_min <- min(rowMeans(beta_temp))
    b_max <- max(rowMeans(beta_temp))
    #09.400
    b_09_400_median <- mean(rowMedians(beta_temp_C09))
    b_09_400_mean <- mean(rowMeans(beta_temp_C09))
    b_09_400_mean_sd <- mean(rowSds(beta_temp_C09))
    b_09_400_sd_mean <- sd(rowMeans(beta_temp_C09))
    b_09_400_mean_CV <-  mean(rowSds(beta_temp_C09)/rowMeans(beta_temp_C09))
    b_09_400_CV_mean <- sd(rowMeans(beta_temp_C09))/mean(rowMeans(beta_temp_C09))
    b_09_400_min <- min(rowMeans(beta_temp_C09))
    b_09_400_max <- max(rowMeans(beta_temp_C09))
    #09.2800
    b_80_400_median <- mean(rowMedians(beta_temp_C80))
    b_80_400_mean <- mean(rowMeans(beta_temp_C80))
    b_80_400_mean_sd <- mean(rowSds(beta_temp_C80))
    b_80_400_sd_mean <- sd(rowMeans(beta_temp_C80))
    b_80_400_mean_CV <-  mean(rowSds(beta_temp_C80)/rowMeans(beta_temp_C80))
    b_80_400_CV_mean <- sd(rowMeans(beta_temp_C80))/mean(rowMeans(beta_temp_C80))
    b_80_400_min <- min(rowMeans(beta_temp_C80))
    b_80_400_max <- max(rowMeans(beta_temp_C80))
    #80.400
    b_09_2800_median <- mean(rowMedians(beta_temp_E09))
    b_09_2800_mean <- mean(rowMeans(beta_temp_E09))
    b_09_2800_mean_sd <- mean(rowSds(beta_temp_E09))
    b_09_2800_sd_mean <- sd(rowMeans(beta_temp_E09))
    b_09_2800_mean_CV <-  mean(rowSds(beta_temp_E09)/rowMeans(beta_temp_E09))
    b_09_2800_CV_mean <- sd(rowMeans(beta_temp_E09))/mean(rowMeans(beta_temp_E09))
    b_09_2800_min <- min(rowMeans(beta_temp_E09))
    b_09_2800_max <- max(rowMeans(beta_temp_E09))
    #80.2800
    b_80_2800_median <- mean(rowMedians(beta_temp_E80))
    b_80_2800_mean <- mean(rowMeans(beta_temp_E80))
    b_80_2800_mean_sd <- mean(rowSds(beta_temp_E80))
    b_80_2800_sd_mean <- sd(rowMeans(beta_temp_E80))
    b_80_2800_mean_CV <-  mean(rowSds(beta_temp_E80)/rowMeans(beta_temp_E80))
    b_80_2800_CV_mean <- sd(rowMeans(beta_temp_E80))/mean(rowMeans(beta_temp_E80))
    b_80_2800_min <- min(rowMeans(beta_temp_E80))
    b_80_2800_max <- max(rowMeans(beta_temp_E80))
    if(i == 1){
      DNAm_gene_summary <- data.frame(gene_id=meta_temp$gene_id[1],size=as.numeric(size),N=nrow(beta_temp),CpG_Range,Span,
                             b_median,b_mean,b_mean_sd,b_sd_mean,b_mean_CV,b_CV_mean,b_min,b_max,
                             b_09_400_median,b_09_400_mean,b_09_400_mean_sd,b_09_400_sd_mean,b_09_400_mean_CV,b_09_400_CV_mean,b_09_400_min,b_09_400_max,
                             b_09_2800_median,b_09_2800_mean,b_09_2800_mean_sd,b_09_2800_sd_mean,b_09_2800_mean_CV,b_09_2800_CV_mean,b_09_2800_min,b_09_2800_max,
                             b_80_400_median,b_80_400_mean,b_80_400_mean_sd,b_80_400_sd_mean,b_80_400_mean_CV,b_80_400_CV_mean,b_80_400_min,b_80_400_max,
                             b_80_2800_median,b_80_2800_mean,b_80_2800_mean_sd,b_80_2800_sd_mean,b_80_2800_mean_CV,b_80_2800_CV_mean,b_80_2800_min,b_80_2800_max,
                             stringsAsFactors = FALSE)
    }else{
      DNAm_gene_summary <- DNAm_gene_summary %>% add_row(gene_id=meta_temp$gene_id[1],size=as.numeric(size),N=nrow(beta_temp),CpG_Range,Span,
                                       b_median,b_mean,b_mean_sd,b_sd_mean,b_mean_CV,b_CV_mean,b_min,b_max,
                                       b_09_400_median,b_09_400_mean,b_09_400_mean_sd,b_09_400_sd_mean,b_09_400_mean_CV,b_09_400_CV_mean,b_09_400_min,b_09_400_max,
                                       b_09_2800_median,b_09_2800_mean,b_09_2800_mean_sd,b_09_2800_sd_mean,b_09_2800_mean_CV,b_09_2800_CV_mean,b_09_2800_min,b_09_2800_max,
                                       b_80_400_median,b_80_400_mean,b_80_400_mean_sd,b_80_400_sd_mean,b_80_400_mean_CV,b_80_400_CV_mean,b_80_400_min,b_80_400_max,
                                       b_80_2800_median,b_80_2800_mean,b_80_2800_mean_sd,b_80_2800_sd_mean,b_80_2800_mean_CV,b_80_2800_CV_mean,b_80_2800_min,b_80_2800_max)
    }
  }
}

# Checking if there is a relationship between number of CpGs and size of gene
plot(DNAm_gene_summary$size~DNAm_gene_summary$N)
# Seems pretty random
# Looking at variation in methylation is related to mean methylation
plot(DNAm_gene_summary$b_sd_mean~DNAm_gene_summary$b_mean)
plot(DNAm_gene_summary$b_mean_CV~DNAm_gene_summary$b_mean)

hist(as.numeric(as.character(DNAm_gene_summary$b_mean)))
# Methylation pattern (which is to say heavily methylated consistent)

locDNAm <- which(!is.na(match(DNAm_gene_summary$gene_id,rownames(rna_cpm))))
dnaM_match <- DNAm_gene_summary[locDNAm,]
dnaM_match <- dnaM_match[order(dnaM_match$gene_id),]
locRNA <- which(!is.na(match(rownames(rna_cpm),DNAm_gene_summary$gene_id,)))
rna_red <- rna_cpm[locRNA,]
i<-1
for(i in 1:nrow(rna_red)){
  print(paste0("Processing gene...",i," of ",nrow(rna_red)))
    # All samples
    rna_median <- median(rna_red[i,])
    rna_mean <- mean(rna_red[i,])
    rna_sd <- sd(rna_red[i,])
    rna_CV <-  sd(rna_red[i,])/mean(rna_red[i,])
    rna_min <- min(rna_red[i,])
    rna_max <- max(rna_red[i,])
    # 09.400
    rna_C09_median <- median(rna_red[i,meta_sample$SFV == "09.400"])
    rna_C09_mean <- mean(rna_red[i,meta_sample$SFV == "09.400"])
    rna_C09_sd <- sd(rna_red[i,meta_sample$SFV == "09.400"])
    rna_C09_CV <-  sd(rna_red[i,meta_sample$SFV == "09.400"])/mean(rna_red[i,meta_sample$SFV == "09.400"])
    rna_C09_min <- min(rna_red[i,meta_sample$SFV == "09.400"])
    rna_C09_max <- max(rna_red[i,meta_sample$SFV == "09.400"])
    # 80.400
    rna_C80_median <- median(rna_red[i,meta_sample$SFV == "80.400"])
    rna_C80_mean <- mean(rna_red[i,meta_sample$SFV == "80.400"])
    rna_C80_sd <- sd(rna_red[i,meta_sample$SFV == "80.400"])
    rna_C80_CV <-  sd(rna_red[i,meta_sample$SFV == "80.400"])/mean(rna_red[i,meta_sample$SFV == "80.400"])
    rna_C80_min <- min(rna_red[i,meta_sample$SFV == "80.400"])
    rna_C80_max <- max(rna_red[i,meta_sample$SFV == "80.400"])
    # 09.2800
    rna_E09_median <- median(rna_red[i,meta_sample$SFV == "09.2800"])
    rna_E09_mean <- mean(rna_red[i,meta_sample$SFV == "09.2800"])
    rna_E09_sd <- sd(rna_red[i,meta_sample$SFV == "09.2800"])
    rna_E09_CV <-  sd(rna_red[i,meta_sample$SFV == "09.2800"])/mean(rna_red[i,meta_sample$SFV == "09.2800"])
    rna_E09_min <- min(rna_red[i,meta_sample$SFV == "09.2800"])
    rna_E09_max <- max(rna_red[i,meta_sample$SFV == "09.2800"])
    # 80.2800
    rna_E80_median <- median(rna_red[i,meta_sample$SFV == "80.2800"])
    rna_E80_mean <- mean(rna_red[i,meta_sample$SFV == "80.2800"])
    rna_E80_sd <- sd(rna_red[i,meta_sample$SFV == "80.2800"])
    rna_E80_CV <-  sd(rna_red[i,meta_sample$SFV == "80.2800"])/mean(rna_red[i,meta_sample$SFV == "80.2800"])
    rna_E80_min <- min(rna_red[i,meta_sample$SFV == "80.2800"])
    rna_E80_max <- max(rna_red[i,meta_sample$SFV == "80.2800"])
    if(i == 1){
      rna_gene_summary <- data.frame(rna_median,rna_mean,rna_sd,rna_CV,rna_min,rna_max,
                            rna_C09_median,rna_C09_mean,rna_C09_sd,rna_C09_CV,rna_C09_min,rna_C09_max,
                            rna_C80_median,rna_C80_mean,rna_C80_sd,rna_C80_CV,rna_C80_min,rna_C80_max,
                            rna_E09_median,rna_E09_mean,rna_E09_sd,rna_E09_CV,rna_E09_min,rna_E09_max,
                            rna_E80_median,rna_E80_mean,rna_E80_sd,rna_E80_CV,rna_E80_min,rna_E80_max,
                            stringsAsFactors = FALSE)
    }else{
      rna_gene_summary <- rna_gene_summary %>% add_row(rna_median,rna_mean,rna_sd,rna_CV,rna_min,rna_max,
                                     rna_C09_median,rna_C09_mean,rna_C09_sd,rna_C09_CV,rna_C09_min,rna_C09_max,
                                     rna_C80_median,rna_C80_mean,rna_C80_sd,rna_C80_CV,rna_C80_min,rna_C80_max,
                                     rna_E09_median,rna_E09_mean,rna_E09_sd,rna_E09_CV,rna_E09_min,rna_E09_max,
                                     rna_E80_median,rna_E80_mean,rna_E80_sd,rna_E80_CV,rna_E80_min,rna_E80_max)
    }
}
## RNA CV vs RNA Expresions
plot(1/rna_gene_summary$rna_CV~log2(rna_gene_summary$rna_mean))
plot(1/rna_gene_summary$rna_C09_CV~log2(rna_gene_summary$rna_C09_mean),col="lightblue")
points(1/rna_gene_summary$rna_C80_CV~log2(rna_gene_summary$rna_E09_mean),col="tomato")

ggplot(rna_gene_summary,aes(y=1/rna_CV,x=log2(rna_mean),colour=as.numeric(dnaM_match$b_mean))) + geom_point()

## Mean DNA methylation vs. mean gene expression
# All Samples
plot(log2(rna_gene_summary$rna_mean)~dnaM_match$b_mean)
ggplot(rna_gene_summary,aes(y=log2(rna_mean),x=dnaM_match$b_mean)) + geom_point() + geom_smooth()
ggplot(rna_gene_summary,aes(y=log2(rna_mean),x=dnaM_match$b_mean)) + geom_point() + geom_smooth() + xlim(0.5,1)
# Control : 09
plot(log2(rna_gene_summary$rna_C09_mean)~dnaM_match$b_09_400_mean,col="lightblue")
# Exposed : 09
points(log2(rna_gene_summary$rna_E09_mean)~dnaM_match$b_09_2800_mean,col="tomato")
# Diff Control vs Exposed: 09
plot(c(log2(rna_gene_summary$rna_C09_mean)-log2(rna_gene_summary$rna_E09_mean))~
     c(dnaM_match$b_09_400_mean-dnaM_match$b_09_2800_mean))
## Mean DNA Methylation vs. CV gene expression
lm_mean_CV <- lm(rna_gene_summary$rna_CV~dnaM_match$b_mean)
summary(lm_mean_CV)
plot(rna_gene_summary$rna_CV~dnaM_match$b_mean)
abline(lm_mean_CV,col="red",lwd=2)
# by treatment : 09
plot(rna_gene_summary$rna_C09_CV~dnaM_match$b_09_400_mean,col="lightblue")
points(rna_gene_summary$rna_E09_CV~dnaM_match$b_09_2800_mean,col="tomato")

## Plot RNA expression difference between C and E at day 9 compared to change in DNAm
diff_exp <- rna_gene_summary$rna_E09_mean-rna_gene_summary$rna_C09_mean
diff_CV <- rna_gene_summary$rna_E09_CV-rna_gene_summary$rna_C09_CV
diff_beta <- dnaM_match$b_09_2800_mean-dnaM_match$b_09_400_mean 
plot(log2(diff_exp)~diff_beta)
plot(diff_CV~diff_beta)

# Checking patterns for genes with > 10 CpGs
dnam_temp <- dnaM_match[dnaM_match$N > 10,]
rna_temp <- rna_gene_summary[dnaM_match$N > 10,]
# Control : 09
plot(log2(rna_temp$rna_C09_mean)~dnam_temp$b_09_400_mean,col="lightblue")
# Exposed : 09
points(log2(rna_temp$rna_E09_mean)~dnam_temp$b_09_2800_mean,col="tomato")

# Control : 09
plot(rna_temp$rna_C09_CV~dnam_temp$b_09_400_mean,col="lightblue")
# Exposed : 09
points(rna_temp$rna_E09_CV~dnam_temp$b_09_2800_mean,col=alpha("tomato",0.4))




