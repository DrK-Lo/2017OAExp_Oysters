Differential Expression with Salmon Gene level outputs with DESeq2
================

### **Data**

``` r
# Meta Data
model<-read.csv("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/metadata_cvirginica_rna_meta.csv", header=TRUE)
model$Treatment <- as.factor(model$treatment)
model$Time <- as.character(model$timepoint)
model$Time[model$Time == "3"] <- "09"
model$Time[model$Time == "6"] <- "80"
model$Time <- as.factor(model$Time)
model$Pop <- as.factor(model$population)
model$Lane <- as.factor(model$lane)
model$SFV <-interaction(model$Time,model$Treatment) # Creates single factor variable for combination of time and treatment

# DESeq Data Object created directly from tximport using `DESeqDataSetFromTximport()` function.
g_obj <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/salmon_pipeline/run20190610_counts/geneDESeqDataObj.RData")

# Additional info about each gene
trans <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/transcriptome_table.RData")
```

**Outdated Count Matrix approach**

``` r
### Gene-level abundances - scaled by isoform length
# Outdated matrix
ga <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/salmon_pipeline/run20190610_counts/geneMatrixAbundance_default.RData")
# # Round counts before create DESeq matrix
ga_round <- round(ga,digits=0)
# Removes just the smallest stuff to make the DESeq analysis run faster (this is not necessary as DESeq also has a filter step built into the result() function that will also remove low count genes from the results)  
ga_filter <- ga_round[rowSums(ga_round)>1,]
dim(ga_filter)
```

### **Full Model** ~ Lane + Pop + shelf + Treatment + Time + Treatment:Time

``` r
# Full Model
# full <- DESeqDataSetFromMatrix(countData = ga_round,
#                               colData = model,
#                               design = ~ Lane + Pop + Treatment + Time + Treatment:Time)

#Perform the DESeq analysis GLM with negative bionomial model and wald test
design(g_obj) <- ~ Lane + Pop + shelf + Treatment + Time + Treatment:Time
g_full_Wald <- DESeq(g_obj,minmu = 80,test = "Wald")
# Interaction
g_full_inter <- results(g_full_Wald)
summary(g_full_inter)
g_pVal_interaction <- g_full_inter$padj
# Time
g_full_Time <- results(g_full_Wald,contrast=c("Time","09","80"))
summary(g_full_Time)
g_pVal_time <- g_full_Time$padj
# Treatment
g_full_Treatment <- results(g_full_Wald,contrast=c("Treatment","400","2800"))
summary(g_full_Treatment)
g_pVal_treatment <- g_full_Treatment$padj
```

#### Combining the summary statistics of expression and DE (log2Fold change + adjP) for each contrast.

``` r
pvals_df <- data.frame(location=row.names(ga),
                         baseMean_Expression = res_full$baseMean,
                         Interaction_Log2_FoldChange=res_full$log2FoldChange,
                         Time_Log2_FoldChange=res_full_Time$log2FoldChange,
                         Treatment_Log2_FoldChange=res_full_Treatment$log2FoldChange,
                         Interaction_P=res_full$padj,
                         Time_P=res_full_Time$padj,
                         Treatment_P=res_full_Treatment$padj)

gene_GE <- merge(trans,pvals_df,by="location")
gene_GE <- gene_GE[!duplicated(gene_GE$location), ]
saveRDS(gene_GE,"/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/salmon_pipeline/run20190610/gene_full_DESeq2Results.RData")
```

**Top Twenty Genes Associated with Treatment**

``` r
gene_GE_trt <- gene_GE[order(gene_GE$Treatment_P),]
kable(gene_GE_trt[1:20,]) %>%
  kable_styling()
```

**Top Twenty Genes Associated with Time**

``` r
gene_GE_t <- gene_GE[order(gene_GE$Time_P),]
kable(gene_GE_t[1:20,]) %>%
  kable_styling()
```

**Top Twenty Genes Associated with Time:Treatment**

``` r
gene_GE_ttrt <- gene_GE[order(gene_GE$Interaction_P),]
kable(gene_GE_ttrt[1:20,]) %>%
  kable_styling()
```

### **Alt Full Model** ~ Lane + Pop + Shelf + SFV

  - SFV = a groups factor with each Time+Treatment as a separate level.
    This allows specific contrasts for each time point and treatment.

<!-- end list -->

``` r
design(g_obj) <- ~ Lane + Pop + shelf + SFV
g_SFV_Wald <- DESeq(g_obj,minmu = 80,test = "Wald")

# Timepoint 9
g_SFV_TP1 <- results(g_SFV_Wald,contrast=c("SFV","09.400","09.2800"))
summary(g_SFV_TP1)
g_SFV_pVal_TP1 <- g_SFV_TP1$padj
# Timepoint 80
g_SFV_TP2 <- results(g_SFV_Wald,contrast=c("SFV","80.400","80.2800"))
summary(g_SFV_TP2)
g_SFV_pVal_TP2 <- g_SFV_TP2$padj
# Tank affect (looking at both controls)
g_SFV_C <- results(g_SFV_Wald,contrast=c("SFV","09.400","80.400"))
summary(g_SFV_C)
g_SFV_pVal_C <- g_SFV_C$padj
# Tank affect (duration effect in OA)
g_SFV_E <- results(g_SFV_Wald,contrast=c("SFV","09.2800","80.2800"))
summary(g_SFV_E)
g_SFV_pVal_E <- g_SFV_E$padj
```

#### Combining the summary statistics of expression and DE (log2Fold change + adjP) for each contrast.

``` r
pvals_df <- data.frame(location=row.names(ga),
                         baseMean_Expression = g_SFV_TP1$baseMean,
                         TP1_Log2_FoldChange=g_SFV_TP1$log2FoldChange,
                         TP2_Log2_FoldChange=g_SFV_TP2$log2FoldChange,
                         Diff_Log2_FoldChange=c(g_SFV_TP1$log2FoldChange-g_SFV_TP2$log2FoldChange),
                         Baseline_Log2_FoldChange=g_SFV_C$log2FoldChange,
                         TP1_P=g_SFV_pVal_TP1,
                         TP2_P=g_SFV_pVal_TP2,
                         Duration_Baseline = g_SFV_pVal_C)

gene_GE <- merge(trans,pvals_df,by="location")
# Removed redudant rows (the transcript variant rows for each gene since this is gene level analysis)
gene_GE <- gene_GE[!duplicated(gene_GE$location), ]
#Remove unwanted columns
gene_GE <- gene_GE[-c(3,5)]
# Rank data frame based on lowest pvalues for either TP1 or TP2
P_val <- cbind(gene_GE$TP1_P,gene_GE$TP2_P) 
P_val[is.na(P_val)] <- 1
P_val_rank <- apply(P_val,1,min)
# Distribution of min P_val (between two timepoints) with p_vals >0.95 removed (most of them)
hist(P_val_rank[P_val_rank<.95])
gene_GE_order <- gene_GE[order(P_val_rank),]
write.csv(gene_GE_top50 <- gene_GE_order[1:50,],"/home/downeyam/Github/2017OAExp_Oysters/results/Salmon_Gene/run20190610_geneList_DESeq2Results_top50.csv")
saveRDS(gene_GE_order,"/home/downeyam/Github/2017OAExp_Oysters/results/Salmon_Gene/run20190610_geneList_DESeq2Results.RData")
```

#### Kable Table Code for top 50 DEGs from either timepoint

(work in
progress)

``` r
labels <- c("Location","Gene ID","Predicted Function","mean Expression","Day 09 : Treatment FoldChange (log2)","Day 80 : Treatment FoldChange (log2)","Baseline FoldChange (log2)","Day 80 P","Day 09 P","Baseline P")
names(gene_GE_order) <- labels

tb <- kable(gene_GE_order[1:50,],format = "latex") %>%
  kable_styling
kableExtra::save_kable(tb,file="/home/downeyam/Github/2017OAExp_Oysters/notebook/img/output2.html")
kableExtra::save_kable(tb,file="/home/downeyam/Github/2017OAExp_Oysters/notebook/img/output2.pdf")
```

### Timepoint Specific DE analysis

**Description**: Here we look are if there is differential expression
due to treatment for each timepoint separately (this is perhaps not the
best way to handle the data, but appears to be the most common approach
for timeseries RNAseq datasets). Also we can look at DE using two
different approaches. The `test = "Wald"` test, which is the default
above and is based on using a GLM with a logarithmic link (negative
binomial) to test for a significant effect of the variable of interest.
Alternatively, I also ran the `test = "LRT"` test, which performs a
likelihood ratio test where the two models being tested in this case are
either with (full model) or without (simple model) treatment (PCO2).

### DEGs specifically for individuals from **TP9**

``` r
ga_tp9 <-  ga_filter[,model$Time == "09"]
dim(ga_tp9)
head(ga_tp9)
model_tp9 <- model[model$Time == "09",]
dim(model_tp9)
tp9 <- DESeqDataSetFromMatrix(countData = ga_tp9,
                              colData = model_tp9,
                              design = ~ Lane + Pop + Treatment)

# Wald test using negative binomial
dds_tp9_Wald <- DESeq(tp9,minmu = 80,test = "Wald")
# Likelihood ratio test
dds_tp9_LTR <- DESeq(tp9,minmu=80,test="LRT",reduced = ~ Lane + Pop)

# Treatment at timepoint 9 with Wald
res_tp9 <- results(dds_tp9_Wald,contrast=c("Treatment","400","2800"))
summary(res_tp9)
pVal_tp9 <- res_tp9$padj

# Treatment at timepoint 9 with LRT
res_tp9_lrt <- results(dds_tp9_LTR)
summary(res_tp9_lrt)
pVal_tp9_lrt <- res_tp9_lrt$padj

pvals_tp9_df <- data.frame(location=row.names(ga_filter),wald_p=pVal_tp9,lrt_p=pVal_tp9_lrt)

gene_GE_tp9 <- merge(trans,pvals_tp9_df,by="location")
gene_GE_tp9 <- gene_GE_tp9[order(gene_GE_tp9$wald_p),]
kable(gene_GE_tp9[1:50,]) %>%
  kable_styling()
```

### DEGs specifically for individuals from **TP80**

``` r
ga_tp80 <-  ga_filter[,model$Time == "80"]
dim(ga_tp80)
head(ga_tp80)
model_tp80 <- model[model$Time == "80",]
dim(model_tp80)
tp80 <- DESeqDataSetFromMatrix(countData = ga_tp80,
                              colData = model_tp80,
                              design = ~ Lane + Pop + Treatment)

# Wald test using negative binomial
dds_tp80_Wald <- DESeq(tp80,minmu = 80,test = "Wald")
# Likelihood ratio test
dds_tp80_LTR <- DESeq(tp80,minmu=80,test="LRT",reduced = ~Lane + Pop)

# Treatment at timepoint 9 with Wald
res_tp80 <- results(dds_tp80_Wald,contrast=c("Treatment","400","2800"))
summary(res_tp80)
pVal_interaction <- res_tp80$padj

# Treatment at timepoint 9 with LRT
res_tp80_lrt <- results(dds_tp80_LTR)
summary(res_tp80_lrt)
pVal_interaction_lrt <- res_tp80_lrt$padj

pvals_tp80_df <- data.frame(location=row.names(gc),wald_p=pVal_tp9,lrt_p=pVal_tp9_lrt)

gene_GE_tp80 <- merge(trans,pvals_tp80_df,by="location")
gene_GE_tp80 <- gene_GE_tp80[order(gene_GE_tp80$wald_p),]
kable(gene_GE_tp80[1:50,]) %>%
  kable_styling()
```
