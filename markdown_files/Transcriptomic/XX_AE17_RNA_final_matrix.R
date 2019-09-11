# Libraries
library(edgeR)
library(limma)
library(DESeq2)
library(dplyr)
library(variancePartition)

#### READ IN DATA ######

#### STAR counts #####
gcounts <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/STAR_mapping/STAR_quantMode_output/STAR_count/readCountMatrix.Rdata")
#### Salmon counts ####
# Counts (not in TPM) : ls = length scaled
s_ls_counts <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/STAR_mapping/salmon_output/salmon_length_ScaledTPM_geneMatrixCounts.RData")
s_counts <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/STAR_mapping/salmon_output/salmon_scaledTPM_geneMatrixCounts.RData")
# Abundance estimates after transformed into TPM (i.e. each sample counts sum to 1 million)
s_ls_abundance <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/STAR_mapping/salmon_output/salmon_length_ScaledTPM_geneMatrixAbundance.RData")
s_abundance <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/STAR_mapping/salmon_output/salmon_scaledTPM_geneMatrixAbundance.RData")
#### RSEM counts ####
RSEM <-  readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/STAR_mapping/RSEM_output/RSEM_gene_Summary.Rdata")
# Separate out RSEM counts and rename rows with LOC ID
rsem_c <- RSEM$Count
#### Original G Matrix ####
# Original count matrix from STAR with HT-seq
gMat <- read.delim("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/C_virginica_gene_count_final.txt",sep = " ")
# Reduce to only include overlapping genes 
gMatrv <- gMat[na.exclude(match(rownames(s_counts),rownames(gMat))),]
na.

#### Transcript File ####
tran <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/STAR_gnomon_tximportGeneFile.RData")
# Gene File
gene <- tran[!duplicated(tran$GENEID),]

  #### Meta Data ####
  meta <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/meta/metadata_20190718.RData")
  meta$sampl_nameSimple <- substr(meta$sample_name,start = 4,stop=9)
  #Creat new factor levels (one for each level combination)
  meta$SFVrn <- as.factor(paste0("D",meta$SFV))
  meta$Sample_Index <- as.factor(meta$sample_index)
  meta$TankID <- as.factor(meta$tankID)
  
#### Data Manipulation ####
# Separate out STAR and RSEM counts and rename rows with LOC ID
loc_star <- gene$GENEID[match(row.names(gcounts),gene$gene_id)]
row.names(gcounts) <- loc_star
loc_rsem <- gene$GENEID[match(row.names(rsem_c),gene$gene_id)]
row.names(rsem_c) <- loc_rsem

# Label columns for salmon based matrices
colnames(s_ls_counts) <- colnames(rsem_c)
colnames(s_counts) <- colnames(rsem_c)

# Reorder different matrixes to perform simpl correlation
star_counts_order <- gcounts[order(rownames(gcounts)),]
s_ls_counts_order <- s_ls_counts[order(rownames(s_ls_counts)),]
s_counts_order <- s_counts[order(rownames(s_counts)),]
rsem_order <- rsem_c[order(rownames(rsem_c)),]

#### Basic Comparisons between mapping and quantification approaches ####

## Genes missing from the current count matrices compared to original matrix by Brett
(na_gMat <- setdiff(rownames(gMat),rownames(gMatrv)))
kableExtra::kable(gMat_missing <- gMat[na.exclude(na_gMat,match(rownames(gMat))),]) %>% kableExtra::kable_styling()

## Summary of counts (STAR v Salmon length scaled v Salmon v RSEM)
star_counts_sum <- colSums(star_counts_order)
s_ls_counts_sum <- colSums(s_ls_counts_order)
s_counts_sum <- colSums(s_counts_order)
rsem_sum <- colSums(rsem_order)
cNames <- c("STAR_quantMode","Salmon_lengthScaled","Salmon","RSEM")
sum_mat <- cbind(star_counts_sum,s_ls_counts_sum,s_counts_sum,rsem_sum)
colnames(sum_mat) <- cNames
rownames(sum_mat) <- colnames(star_counts_order)
final_sum_mat <- cbind(apply(sum_mat,2,mean),apply(sum_mat,2,sd),
                       apply(sum_mat,2,min),apply(sum_mat,2,max))
rownames(final_sum_mat) <- cNames
colnames(final_sum_mat) <- c("mean","sd","min","max")
kableExtra::kable(final_sum_mat) %>% kableExtra::kable_styling()


## Pearson correlation between methods (STAR v Salmon length scaled v Salmon v RSEM)
starVs <- NULL
starVsls <- NULL
starVrsem <- NULL
sVsls <- NULL
sVrsem <- NULL
slsVrsem <- NULL

for(i in 1:24){
  starVs[i] <- cor.test(star_counts_order[,i],s_counts_order[,i])$estimate
  starVsls[i] <- cor.test(s_ls_counts_order[,i],star_counts_order[,i])$estimate
  starVrsem[i] <- cor.test(star_counts_order[,i],rsem_order[,i])$estimate
  sVsls[i] <- cor.test(s_ls_counts_order[,i],s_counts_order[,i])$estimate
  slsVrsem[i] <- cor.test(s_ls_counts_order[,i],rsem_order[,i])$estimate
  sVrsem[i] <- cor.test(s_counts_order[,i],rsem_order[,i])$estimate
}
cor_mat <- rbind(starVs,starVsls,starVrsem,sVsls,sVrsem,slsVrsem)
corr_mean <- cbind(apply(cor_mat,1,mean),apply(cor_mat,1,sd))
colnames(corr_mean) <- c("mean","sd")
kableExtra::kable(corr_mean) %>% kableExtra::kable_styling()
###### LIMMA - VOOM PIPELINE #######

#### Filtering and normalization - limma - all samples - RSEM/salmon_lengthscaled/salmon ####
# Create DGEList
g <- DGEList(gcounts,genes = gene) # Star quantMode
s_ls <- DGEList(s_ls_counts_order,genes = gene) # length scaled counts - salmon
s <- DGEList(s_counts_order,genes = gene) # counts - salmon
rsem <- DGEList(rsem_order,genes = gene) # counts - rsem

# Identify positions to filter reads
# Removes genes that have a expression <1count per million across all individuals
# 24 equal the number of samples (so at least half need to have a count)
keep_g <- rowSums(cpm(g)>1) >= (0.5 * 24)
keep_s_ls <- rowSums(cpm(s_ls)>1) >= (0.5 * 24)
keep_s <- rowSums(cpm(s)>1) >= (0.5 * 24)
keep_rsem <- rowSums(cpm(rsem)>1) >= (0.5 * 24)
# Outdated approach from edgeR for filtering uninformative genes (those with little expression)
#keep_s_ls <- filterByExpr(s_ls)
#keep_s <- filterByExpr(s)
#keep_rsem <- filterByExpr(rsem)

# Filter reads
g_red <- g[keep_g, ]
s_ls_red <- s_ls[keep_s_ls, ]
s_red <- s[keep_s, ]
rsem_red <- rsem[keep_rsem, ]
# Check number of genes lost
dim(g_red)
dim(s_ls_red)
dim(s_red)
dim(rsem_red)

# Calculate normalization factor
g_norm <- calcNormFactors(g_red)
s_ls_norm <- calcNormFactors(s_ls_red)
s_norm <- calcNormFactors(s_red)
rsem_norm <- calcNormFactors(rsem_red)

# Basic design matrix to start
design <- model.matrix(~0 + Treatment*Time, data = meta)
#Performing voom function with MDS plots
g_voom <- voom(g_norm,design,plot = TRUE)
plotMDS(g_voom,col=as.numeric(meta$Treatment),main="Treatment : STAR Quantmode")
s_ls_voom <- voom(s_ls_norm, design,plot = TRUE)
plotMDS(s_ls_voom,col=as.numeric(meta$Treatment),main="Treatment : salmon length scaled")
s_voom <- voom(s_norm, design,plot=TRUE)
plotMDS(s_voom,col=as.numeric(meta$Treatment),main="Treatment : salmon")
rsem_voom <- voom(rsem_norm, design,plot=TRUE)
plotMDS(rsem_voom,col=as.numeric(meta$Treatment),main="Treatment : RSEM")
plotMDS(rsem_voom,col=as.numeric(meta$Time),main="Time : RSEM")
plotMDS(rsem_voom,col=as.numeric(meta$Pop),main="Population : RSEM")

#### Basic DE Analysis (no random facs or contrasts) - limma - all samples  - RSEM/salmon_lengthscaled/salmon ####
## Top DEGs for each sample 
#star quantMode
g_fit <- lmFit(g_voom, design)
g_fit <- eBayes(g_fit)
g_top <- topTable(g_fit, coef=ncol(design))
kableExtra::kable(g_top) %>% kableExtra::kable_styling()
#salmon lengthscaled
s_ls_fit <- lmFit(s_ls_voom, design)
s_ls_fit <- eBayes(s_ls_fit)
s_ls_top <- topTable(s_ls_fit, coef=ncol(design))
kableExtra::kable(s_ls_top) %>% kableExtra::kable_styling()
#salmon
s_fit <- lmFit(s_voom, design)
s_fit <- eBayes(s_fit)
s_top <- topTable(s_fit, coef=ncol(design))
kableExtra::kable(s_top) %>% kableExtra::kable_styling()
#rsem
rsem_fit <- lmFit(rsem_voom, design)
rsem_fit <- eBayes(rsem_fit)
rsem_top <- topTable(s_fit, coef=ncol(design))
kableExtra::kable(rsem_top) %>% kableExtra::kable_styling()
## Really not significant so far....

#### Basic DE Analysis (Timepoints separate (09 and 80)) - limma - RSEM ####

# Day 09
meta_09 <- meta[meta$Day == 9,]
design <- model.matrix(~Treatment, data = meta_09)

sp_rsem <- match(meta_09$sampl_nameSimple,colnames(rsem_norm))
rsem_norm_09 <- rsem_norm[,sp_rsem]
rsem_norm_09 <- calcNormFactors(rsem_norm_09)

rsem_09_voom <- voom(rsem_norm_09, design,plot = TRUE)
plotMDS(rsem_09_voom,col=as.numeric(meta_09$Treatment))

rsem_09_fit <- lmFit(rsem_09_voom,design)
rsem_09_fit <- eBayes(rsem_09_fit)
rsem_09_top <- topTable(rsem_09_fit, coef=ncol(design),sort.by = "p")
kableExtra::kable(rsem_09_top) %>% kableExtra::kable_styling()

# Day 80
meta_80 <- meta[meta$Day == 80,]
design <- model.matrix(~Treatment, data = meta_80)

sp_rsem <- match(meta_80$sampl_nameSimple,colnames(rsem_norm))
rsem_norm_80 <- rsem_norm[,sp_rsem]
rsem_norm_80 <- calcNormFactors(rsem_norm_80)

rsem_80_voom <- voom(rsem_norm_80, design,plot = TRUE)
plotMDS(rsem_80_voom,col=as.numeric(meta_80$Treatment))

rsem_80_fit <- lmFit(rsem_80_voom, design)
rsem_80_fit <- eBayes(rsem_80_fit)
rsem_80_top <- topTable(rsem_80_fit, coef=ncol(design))
kableExtra::kable(rsem_80_top) %>% kableExtra::kable_styling()

# Removing sample 17005 from day 80 analysis as potential outlier
# meta_80 <- meta_80[meta_80$sampl_nameSimple != "17005",]
# design <- model.matrix(~Treatment, data = meta_80)
# 
# sp_rsem <- match(meta_80$sampl_nameSimple,colnames(rsem_red))
# rsem_red_80 <- rsem_red[,sp_rsem]
# rsem_red_80 <- calcNormFactors(rsem_red_80)
# 
# rsem_80_out <- voom(rsem_red_80, design,plot = TRUE)
# plotMDS(rsem_80_out,col=as.numeric(meta_80$Treatment))
# 
# rsem_80_fit <- lmFit(rsem_80_out, design)
# rsem_80_fit <- eBayes(rsem_80_fit)
# topTable(rsem_80_fit, coef=ncol(design),sort.by = "p")

##### DE Analysis (Planned Contrasts) - limma - RSEM ####

#Make new design matrix
design_contr <- model.matrix(~0+SFVrn,data=meta)
#Rename columns
colnames(design_contr) <- levels(meta$SFVrn)
#Create specific contrasts
rsem_contr_mat <- makeContrasts(
  CvE_D9 = D09.2800-D09.400,
  CvE_D80 = D80.2800-D80.400,
  C_D9vD80 = D09.400-D80.400,
  Diff = (D09.2800-D09.400)- (D80.2800-D80.400),
  levels=design_contr
)

# Identify correlation between factors in design contrasts with blocking factor
cor <- duplicateCorrelation(rsem_voom, design_contr, block = meta$tankID)
#Rerun lmFit model w/ blocking factor
rsem_fit_contr <- lmFit(rsem_voom, design_contr,block = meta$tankID,correlation = cor$consensus.correlation)
#Refit with contrasts
rsem_fit2_contr <- contrasts.fit(rsem_fit_contr,rsem_contr_mat)
#Run in eBayes
rsem_fit2_contr <- eBayes(rsem_fit2_contr)
#Examine results
topTable(rsem_fit2_contr)
rsem_results_contr <- decideTests(rsem_fit2_contr,adjust.method = "BH")
vennDiagram(rsem_results_contr)

plotMD(rsem_fit2_contr,column=2)

##### Variance Paritioning (from limma) - RSEM #####

# Checking correlation between variable
form <- ~ TankID + Lane + Treatment + Time + Pop + Sample_Index
C = canCorPairs(form,meta)
plotCorrMatrix( C )

## Handling random effects (batch effects, then performing DE)

modelFit <- fitVarPartModel(rsem_voom,~(1|TankID) + (1|Pop),meta)
rsem_res <-  residuals(modelFit)

var_form <- ~ (1|Treatment) + (1|Time)
rsem_varPartResid <- fitExtractVarPartModel(rsem_res,var_form,meta)
rownames(rsem_varPartResid) <- rownames(rsem_voom)
rsem_sort_varPartResid <- sortCols(sortCols(rsem_varPartResid),)


plotVarPart(rsem_sort_varPartResid,label.angle=60)
head(rsem_sort_varPartResid)
top_time <- rsem_sort_varPartResid[order(rsem_sort_varPartResid$Time,decreasing = TRUE),]
plotPercentBars(top_time[1:10,])
top_trt <- rsem_sort_varPartResid[order(rsem_sort_varPartResid$Treatment,decreasing = TRUE),]
plotPercentBars(top_trt[1:10,])

#### DE Anaysis - DESeq2 - RSEM #####

# Full Mode

g_obj <- DESeqDataSetFromMatrix(countData = round(rsem_order),
                               colData = meta,
                               design = ~ Lane + Pop + Treatment + Time + Treatment:Time)
#Perform the DESeq analysis GLM with negative bionomial model and wald test
#design(g_obj) <- ~ Lane + Pop + shelf + Treatment + Time + Treatment:Time
g_full_Wald <- DESeq(g_obj,minmu = 80,test = "Wald")
g_full_Wald_scaling <- estimateSizeFactors(g_full_Wald)

# Interaction
g_full_inter <- results(g_full_Wald_scaling)
summary(g_full_inter)
g_pVal_interaction <- g_full_inter$padj
# Time
g_full_Time <- results(g_full_Wald_scaling,contrast=c("Time","09","80"))
summary(g_full_Time)
g_pVal_time <- g_full_Time$padj
# Treatment
g_full_Treatment <- results(g_full_Wald_scaling,contrast=c("Treatment","400","2800"))
summary(g_full_Treatment)
g_pVal_treatment <- g_full_Treatment$padj

#### Variance Paritioning (from DESeq2 - RSEM) #####
keep_g_DESeq <- rowSums(fpm(g_full_Wald_scaling)>1) >= (0.5 * 24)
g_DES <- log2(fpm(g_full_Wald_scaling)[keep_g_DESeq,] + 1)

# Single step without first accounting for true random factors
var_form <- ~ (1|Treatment) + (1|Time) +(1|TankID) + (1|Pop)
var_form <- ~ Treatment + Time + 
varPart_DESeq_rsem <- fitExtractVarPartModel(g_DES,var_form,meta)

rownames(varPart_DESeq_rsem) <- rownames(g_DES)
rsem_DESeq_sort_varPartResid <- sortCols(sortCols(varPart_DESeq_rsem),)


plotVarPart(rsem_DESeq_sort_varPartResid,label.angle=60)
