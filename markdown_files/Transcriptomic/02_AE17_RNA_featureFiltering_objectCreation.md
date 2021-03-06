Filtering and Storing Expression Matrices
================
Alan Downey-Wall
10/18/2019

### Purpose

In this step you will filter the feature count matrices (gene,
transcript, and exon) to include only features that only have expression
in at least 5 samples within a treatment-time combination. See
[`AE17_geneFilteringComparison`]() for the full break down on filtering
and normalization criteria. Here we also used the filter matrices to
create DGEList objects which store information about each feature and
can be used by `EdgeR` or `limma` for downstream [differential
analyses]() or for [genome-wide PCA or DAPC approaches]().

### Inputs

  - Raw Count
        Matrices
      - `/input_files/RNA/STAR_mapping/RSEM_output/RSEM_gene_Summary.Rdata`  
      - `/input_files/RNA/STAR_mapping/RSEM_output/RSEM_isoform_Summary.Rdata`
  - Feature Annotation File
      - `/input_files/RNA/references/STAR_gnomon_tximportGeneFile.RData`
  - Sample Meta Data File
      - `/input_files/meta/metadata_20190811.RData`

### Outputs

  - Post-filtering, pre-normalization DGEList
        objects
      - `/results/Transcriptomic/gene_preNormalization_DGEListObj.RData`  
      - `/results/Transcriptomic/transcript_preNormalization_DGEListObj.RData`  
  - Post voom and normalization DGEList
        objects
      - `/results/Transcriptomic/gene_postVoomAndNormalization_DGEListObj.RData`  
      - `/results/Transcriptomic/gene_postVoomAndNormalization_DGEListObj.RData`  
  - Output from differential expression using `limma-voom` linear model
    fit approach
      - `/results/Transcriptomic/gene_EBayesObj.RData`  
      - `/results/transcript_EBayesObj.RData`

**Libraries and setting the path**

``` r
## Packages 
library(matrixStats,quietly = TRUE)
library(edgeR,quietly = TRUE)
library(limma,quietly = TRUE)
library(fdrtool,quietly = TRUE)
library(kableExtra,quietly = TRUE)
library(DESeq2,quietly = TRUE)
library(dplyr,quietly = TRUE)
library(variancePartition,quietly = TRUE)
library(FSA,quietly = TRUE) # to use w arguement in hist()

wd <- "/home/downeyam/Github/2017OAExp_Oysters"
# This should be set to the path for the local version of the `2017OAExp_Oysters` github repo.
```

**Data**  
*Gene and Transcript Expression Matrices*

``` r
#### RSEM counts ####
# Gene matrix
RSEM <-  readRDS(paste0(wd,"/input_files/RNA/STAR_mapping/RSEM_output/RSEM_gene_Summary.Rdata"))
# Separate out RSEM counts and rename rows with LOC ID
rsem_c <- RSEM$Count # Stick with raw gene counts  (not TPM or some other normalizaed count)
rm(RSEM)

# Isoform (transcript) matrix
RSEM_t <-  readRDS(paste0(wd,"/input_files/RNA/STAR_mapping/RSEM_output/RSEM_isoform_Summary.Rdata"))
rsem_t_c <- RSEM_t$Count
rm(RSEM_t)
```

*Transcript annotation file*

``` r
#### Transcript File ####
# Transcript file
tran <- readRDS(paste0(wd,"/input_files/RNA/references/STAR_gnomon_tximportGeneFile.RData"))
# Gene File
gene <- tran[!duplicated(tran$GENEID),]
```

*Sample meta data*

``` r
#### Meta Data ####
meta <- readRDS(paste0(wd,"/input_files/meta/metadata_20190811.RData"))
meta$sampl_nameSimple <- substr(meta$sample_name,start = 4,stop=9)
#Create new factor levels (one for each level combination)
meta$SFVrn <- as.factor(paste0("D",meta$SFV))
meta$Sample_Index <- as.factor(meta$sample_index)
meta$TankID <- as.factor(meta$tankID)
```

*Data Clean-up*

Here I reordered both the count matrices (gene and transcript) and
annotation objects for when I combine them in the DGEList object
downstream.

``` r
#### Data Manipulation ####
# Order genes from annotation list to match the order in the count matrix
gene_order <- gene[order(gene$gene_id),]
identical(rownames(rsem_c),gene_order$gene_id) # TRUE confirms the order matches
```

    ## [1] TRUE

``` r
# Relabel the rows of the count matrix with LOC ID
rownames(rsem_c) <- gene_order$GENEID

# Order transcripts from annotation list to match the order from the count matrix
rsem_t_order <- rsem_t_c[order(rownames(rsem_t_c)),]
tran_order <- tran[order(tran$TXNAME),]
identical(rownames(rsem_t_order),tran_order$TXNAME)
```

    ## [1] TRUE

``` r
# I am not renaming the rows here since the TXNAME will be unique, while the gene LOC will be redundant for isoforms from the same gene. However, the reordering here is still needed when the annotations get added to the DGEList later.

# Reorder by LOCID 
geneC_full <- rsem_c
tranC_full <- rsem_t_order

rm(rsem_c)
rm(rsem_t_order)
```

### Filtering matrices

``` r
# Round to whole counts
geneC_all <- round(geneC_full)
tranC_all <- round(tranC_full)

### Genes ### 
# Breaking down expression coverage by treatment*time combination
#Day 9 Trt 2800
keep_D9.2800 <- rowSums(cpm(geneC_all[,meta$SFVrn=="D09.2800"])>=1) >= 5
sum(keep_D9.2800)
#Day 9 Trt 400
keep_D9.400 <- rowSums(cpm(geneC_all[,meta$SFVrn=="D09.400"])>=1) >= 5
sum(keep_D9.400)
#Day 80 Trt 2800
keep_D80.2800 <- rowSums(cpm(geneC_all[,meta$SFVrn=="D80.2800"])>=1) >= 5
sum(keep_D80.2800)
#Day 80 Trt 400
keep_D80.400 <- rowSums(cpm(geneC_all[,meta$SFVrn=="D80.400"])>=1) >= 5
sum(keep_D80.400)

keep_gene_a2 <- rowSums(cbind(keep_D9.2800,keep_D9.400,
                                 keep_D80.2800,keep_D80.400)) >= 1
paste0(sum(keep_gene_a2)," of ",summary_df$Num_Features[1]," genes kept after filtering (",
       round(sum(keep_gene_a2)/summary_df$Num_Features[1]*100,2),"% remaining).")

# Filter 
geneC_a2 <- geneC_all[keep_gene_a2, ]
gene_final_a2 <- gene_order[keep_gene_a2,]

### Transcripts ###
# Breaking down expression coverage by treatment*time combination
#Day 9 Trt 2800
keep_D9.2800 <- rowSums(cpm(tranC_all[,meta$SFVrn=="D09.2800"])>=1) >= 5
sum(keep_D9.2800)
#Day 9 Trt 400
keep_D9.400 <- rowSums(cpm(tranC_all[,meta$SFVrn=="D09.400"])>=1) >= 5
sum(keep_D9.400)
#Day 80 Trt 2800
keep_D80.2800 <- rowSums(cpm(tranC_all[,meta$SFVrn=="D80.2800"])>=1) >= 5
sum(keep_D80.2800)
#Day 80 Trt 400
keep_D80.400 <- rowSums(cpm(tranC_all[,meta$SFVrn=="D80.400"])>=1) >= 5
sum(keep_D80.400)

keep_tran_a2 <- rowSums(cbind(keep_D9.2800,keep_D9.400,
                                 keep_D80.2800,keep_D80.400)) >= 1
paste0(sum(transcripts_gene_a2)," of ",summary_df$Num_Features[1]," transcripts kept after filtering (",
       round(sum(keep_gene_a2)/summary_df$Num_Features[1]*100,2),"% remaining).")
#Filter 
tranC_a1 <- tranC_all[keep_tran_a2, ]
tran_final <- tran_order[keep_tran_a2, ]

## Create DGEList
dge_gene_a2 <- DGEList(geneC_a2,genes = gene_final_a2) # counts - rsem
dge_tran_a1 <- DGEList(tranC_a1,genes = tran_final) # counts - rsem

# Save initial objects   
saveRDS(dge_gene_a2,paste0(wd,"/results/Transcriptomic/gene_preNormalization_DGEListObj.RData"))
saveRDS(dge_tran_a1,paste0(wd,"/results/Transcriptomic/transcript_preNormalization_DGEListObj.RData"))
```

### Normalization with edgeR

``` r
# Calculate normalization factors for scaling raw lib. size
dge_gene_a2_norm <- calcNormFactors(dge_gene_a2,method = "TMMwsp") # gene - approach 2
dge_tran_a1_norm <- calcNormFactors(dge_tran_a1,method="TMMwsp")   # transcript - approach 1
```

### Create design matrix

  - **design** - Combining time (2 level factor) and treatment (2 level
    factor) into a single 4 level explanatory factor for the purpose of
    making the specific contrasts (planned comparisons easier to code).
    This was a stragey implemented in the limma manual.
      - Single explanatory factory - **SFVrn** - “D09.2800” “D09.400”
        “D80.2800”
“D80.400”

<!-- end list -->

``` r
design <- model.matrix(~0+SFVrn,data=meta) # 0+ is needed here otherwise the first level defaults to 1.
#Rename columns
colnames(design) <- levels(meta$SFVrn)
```

### Contrasts using `design` matrix

Using the **Option 1** matrix it is easier to perform some of the focal
contrasts that I was particularly interested in. Specifically, I was
interested in three primary contrasts + the interaction:  
1\) CvE\_D9 : Difference between ambient and OA at day 9  
2\) CvE\_D80 : Difference between ambient and OA at day 80  
3\) C\_D9D80 : Difference at ambient between day 9 and day 80  
4\) Diff : Interaction between treatment and time

  - In this case the Diff should be equivalent to the interaction
    (`Treatment:Time`) from *Option 2*

<!-- end list -->

``` r
#Create specific contrasts for option 1
contr_mat <- makeContrasts(
  CvE_D9 = D09.2800-D09.400,
  CvE_D80 = D80.2800-D80.400,
  C_D9vD80 = D09.400-D80.400,
  Diff = (D09.2800-D09.400)- (D80.2800-D80.400),
  levels=design
)
```

### Transform and create observational level weights

  - Use `voomWithQualityWeights` from `limma` package
      - Part of limma package. Will transform count data to log2-count
        per million (logCPM), and estimate the mean-variance
        relationship. This will then be leveraged to compute the
        appropriate observation-level weights.  
      - " When a multi-dimensional scaling plot from a designed RNA-seq
        experiment indicates the presenceof outlier samples, it is
        possible to combine the observation-level weighting strategy
        used in voomwith sample-specific quality weights (as described
        in the section above on Array Quality Weights) todown-weight
        outlier samples. This capability is implemented in
        thevoomWithQualityWeightsfunction. "

<!-- end list -->

``` r
## Gene Features 
dge_gene_a2_o1_voom <- voomWithQualityWeights(dge_gene_a2_norm,design,plot = FALSE)
## Transcript Features 
tran_voom <- voomWithQualityWeights(dge_tran_a1_norm, design,plot = FALSE)

# Saving the final transformation of the data (after individual weights)
saveRDS(dge_gene_a2_o1_voom,paste0(wd,"/results/Transcriptomic/gene_postVoomAndNormalization_DGEListObj.RData"))
saveRDS(tran_voom,paste0(wd,"/results/Transcriptomic/tran_postVoomAndNormalization_DGEListObj.RData"))
```

### Identify correlation between factors in design contrasts with blocking factor

  - Use `duplicateCorrelation` from `limma` package
      - Estimate the correlation between technical replicates within the
        experiment. In this case tanks were replicate treatments (6
        tanks per treatment). Each tank has two individuals (one from
        each timepoint). 12 tanks x 2 individuals = 24 samples.

**Note**: Same as `voom` transformation, the 4 gene feature objects are
used and the 1 transcript feature.

``` r
#Genes
gene_a2_o1_cor <- duplicateCorrelation(dge_gene_a2_o1_voom, design, block = meta$tankID)
# Transcript
tran_cor <- duplicateCorrelation(tran_voom, design, block = meta$tankID)

# Ran this model with and without this step. It doesn't seem to impact the outcome all that much.
```

### Run `lmFit` model and `Ebayes` from `limma`

**Fitting Model without contrasts**

``` r
#Genes
lmf_gene_a2_o1 <- lmFit(dge_gene_a2_o1_voom, design,
                        block = meta$tankID,
                        correlation = gene_a2_o1_cor$consensus.correlation)
#Transcripts
lmf_tran <- lmFit(tran_voom, design,
                        block = meta$tankID,
                        correlation = tran_cor$consensus.correlation)
```

**Refitting Option 1 with contrasts**

``` r
# Genes
lmf_gene_a2_o1_wContr <- contrasts.fit(lmf_gene_a2_o1,contr_mat)
#Transcripts 
lmf_tran_wContr <- contrasts.fit(lmf_tran,contr_mat)
```

**Run empiricial bayes protocol**

``` r
#Run in eBayes w contrasts
#Genes
bayes_gene_a2_o1_contr <- eBayes(lmf_gene_a2_o1_wContr,robust=TRUE)
#Transcript
bayes_tran_contr <- eBayes(lmf_tran_wContr,robust=TRUE)

# Saving Files
saveRDS(bayes_gene_a2_o1_contr,paste0(wd,"/results/Transcriptomic/gene_EBayesObj.RData"))
saveRDS(bayes_tran_contr,paste0(wd,"/results/transcript_EBayesObj.RData"))  
# Ready to perform a differential analyses
```
