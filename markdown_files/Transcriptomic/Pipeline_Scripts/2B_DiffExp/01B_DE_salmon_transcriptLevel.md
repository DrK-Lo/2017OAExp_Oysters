Differential Expression with Salmon Transcript level outputs with DESeq2
================

### **Data**

``` r
# Meta Data
model<-read.csv("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/metadata_cvirginica_rna_meta.txt", header=TRUE)
model$Treatment <- as.factor(model$treatment)
model$Time <- as.character(model$timepoint)
model$Time[model$Time == "3"] <- "09"
model$Time[model$Time == "6"] <- "80"
model$Time <- as.factor(model$Time)
model$Pop <- as.factor(model$population)
model$Lane <- as.factor(model$lane)
model$SFV <-interaction(model$Time,model$Treatment) # Creates single factor variable for combination of time and treatment

# Transcript abundances (only needed for timepoint specific analysis)
ta <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/salmon_pipeline/run20190610_counts/tranMatrixCounts.RData")
# Round counts before create DESeq matrix
ta_round <- round(ta,digits=0)
ta_filter <- ta_round[rowSums(ta_round)>80,]

# Additional info about each gene
# DESeq Data Object created directly from tximport using `DESeqDataSetFromTximport()` function.
t_obj <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/salmon_pipeline/run20190610_counts/tranDESeqDataObj.RData")
# Transcriptome annotations
trans <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/transcriptome_table.RData")
```

### **Full Model** - Lane, Pop, Treatment, Time, Treatment:Time as factors

``` r
### Full Model
# Outdated
# full <- DESeqDataSetFromMatrix(countData = ga_round,
#                               colData = model,
#                               design = ~ Lane + Pop + Treatment + Time + Treatment:Time)

t_full_Wald <- DESeq(t_obj,minmu = 80,test = "Wald")

# Interaction
t_res_full <- results(t_full_Wald,independentFiltering = TRUE)
summary(t_res_full)
pVal_interaction <- t_res_full$padj
# Time
t_res_full_Time <- results(t_full_Wald,contrast=c("Time","09","80"))
summary(t_res_full_Time)
t_pVal_time <- t_res_full_Time$padj
# Treatment
t_res_full_Treatment <- results(t_full_Wald,contrast=c("Treatment","400","2800"))
summary(t_res_full_Treatment)
t_pVal_treatment <- t_res_full_Treatment$padj

t_pvals_df <- data.frame(location=row.names(a),
                         baseMean_Expression = t_res_full$baseMean,
                         Interaction_Log2_FoldChange=t_res_full$log2FoldChange,
                         Time_Log2_FoldChange=t_res_full_Time$log2FoldChange,
                         Treatment_Log2_FoldChange=t_res_full_Treatment$log2FoldChange,
                         Interaction_P=t_res_full$padj,
                         Time_P=t_res_full_Time$padj,
                         Treatment_P=t_res_full_Treatment$padj)

tran_GE <- merge(trans,t_pvals_df,by.x="fullID",by.y="location")
saveRDS(tran_GE,"/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/salmon_pipeline/run20180512/transcript_DESeq2Results.RData")
```

**Top Twenty Genes Associated with Treatment**

``` r
tran_GE_trt <- tran_GE[order(tran_GE$Treatment_P),]
kable(tran_GE_trt[1:20,]) %>%
  kable_styling()
```

**Top Twenty Genes Associated with Time**

``` r
tran_GE_t <- tran_GE[order(tran_GE$Time_P),]
kable(tran_GE_t[1:20,]) %>%
  kable_styling()
```

**Top Twenty Genes Associated with Time:Treatment**

``` r
tran_GE_ttrt <- tran_GE[order(tran_GE$Interaction_P),]
kable(tran_GE_ttrt[1:20,]) %>%
  kable_styling()
```

### **Alt Full Model** ~ Lane + Pop + Shelf + SFV

  - SFV = a groups factor with each Time+Treatment as a separate level.
    This allows specific contrasts for each time point and treatment.

<!-- end list -->

``` r
design(t_obj) <- ~ Lane + Pop + shelf + SFV
t_SFV_Wald <- DESeq(t_obj,minmu = 80,test = "Wald")

# Timepoint 9
t_SFV_TP1 <- results(t_SFV_Wald,contrast=c("SFV","09.400","09.2800"))
summary(t_SFV_TP1)
?results
t_SFV_pVal_TP1 <- t_SFV_TP1$padj
# Timepoint 80
t_SFV_TP2 <- results(t_SFV_Wald,contrast=c("SFV","80.400","80.2800"))
summary(t_SFV_TP2)
t_SFV_pVal_TP2 <- t_SFV_TP2$padj
# Tank affect (looking at both controls)
t_SFV_C <- results(t_SFV_Wald,contrast=c("SFV","09.400","80.400"))
summary(t_SFV_C)
t_SFV_pVal_C <- t_SFV_C$padj
# Tank affect (duration effect in OA)
t_SFV_E <- results(t_SFV_Wald,contrast=c("SFV","09.2800","80.2800"))
summary(t_SFV_E)
t_SFV_pVal_E <- t_SFV_E$padj

#Single  Transcript plots
ordered_val_full <- order(t_SFV_TP1$padj)
ordered_FC_full <- t_SFV_TP1$log2FoldChange
out <- t_obj@assays@.xData$`.->data`$counts
# Currently plotting most significant transcript
plotCounts(t_SFV_Wald, gene=ordered_val_full[1], intgroup="SFV",normalized = FALSE)
log2(out[ordered_val_full[1],])
ordered_FC_full[ordered_val_full[1]]
plotCounts(t_SFV_Wald, gene=ordered_val_full[1], intgroup="SFV",normalized = FALSE)
```

#### Combining the summary statistics of expression and DE (log2Fold change + adjP) for each contrast.

``` r
pvals_df <- data.frame(fullID=row.names(ta),
                         baseMean_Expression = t_SFV_TP1$baseMean,
                         TP1_Log2_FoldChange=t_SFV_TP1$log2FoldChange,
                         TP2_Log2_FoldChange=t_SFV_TP2$log2FoldChange,
                         Diff_Log2_FoldChange=c(t_SFV_TP1$log2FoldChange-t_SFV_TP2$log2FoldChange),
                         Baseline_Log2_FoldChange=t_SFV_C$log2FoldChange,
                         TP1_P=t_SFV_pVal_TP1,
                         TP2_P=t_SFV_pVal_TP2,
                         Duration_Baseline = t_SFV_pVal_C)
head(pvals_df)
head(trans)
tran_GE <- merge(trans,pvals_df,by="fullID")
# Removed redudant rows (the transcript variant rows for each gene since this is gene level analysis)
#gene_GE <- gene_GE[!duplicated(gene_GE$location), ]
#Remove unwanted columns
tran_GE <- tran_GE[-c(3,5)]
# Rank data frame based on lowest pvalues for either TP1 or TP2
P_val <- cbind(tran_GE$TP1_P,tran_GE$TP2_P) 
P_val[is.na(P_val)] <- 1
P_val_rank <- apply(P_val,1,min)
# Distribution of min P_val (between two timepoints) with p_vals >0.95 removed (most of them)
hist(P_val_rank[P_val_rank<.95])
tran_GE_order <- tran_GE[order(P_val_rank),]
LOC <- (unique(tran_GE_order[1:75,]$location))
list <- tran_GE_order[,c(9,10)]
list[is.na(list)] <- 1
LOC_order_Min <- apply(list,1,min)
tran_GE_topHits <- tran_GE_order[LOC_order_Min < .105,2]
tran_GE_order$rank <- 1:nrow(tran_GE_order)
tran_sub <- tran_GE_order[!is.na(match(tran_GE_order$location,unique(tran_GE_topHits))),]
tran_final_order <- arrange(tran_sub,location,rank)
write.csv(tran_final_order,"/home/downeyam/Github/2017OAExp_Oysters/results/Salmon_Gene/run20190610_TPCompare_DESeq2Results_topTrans.csv")
saveRDS(tran_final_order,"/home/downeyam/Github/2017OAExp_Oysters/results/Salmon_Gene/run20190610_tranList_DESeq2Results.RData")
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
ta_tp9 <-  ta_filter[,model$Time == "09"]
dim(ta_tp9)
head(ta_tp9)
model_tp9 <- model[model$Time == "09",]
dim(model_tp9)
tp9 <- DESeqDataSetFromMatrix(countData = ta_tp9,
                              colData = model_tp9,
                              design = ~ Lane + Pop + Treatment)

# Wald test using negative binomial
t_tp9_Wald <- DESeq(tp9,minmu = 80,test = "Wald")
# Likelihood ratio test
t_tp9_LTR <- DESeq(tp9,minmu=80,test="LRT",reduced = ~Lane + Pop)

# Treatment at timepoint 9 with Wald
res_tp9 <- results(t_tp9_Wald)
summary(res_tp9)
pVal_tp9 <- res_tp9$padj

# Treatment at timepoint 9 with LRT
res_tp9_lrt <- results(t_tp9_LTR)
summary(res_tp9_lrt)
pVal_tp9_lrt <- res_tp9_lrt$padj

pvals_tp9_df <- data.frame(location=row.names(ta_filter),wald_p=pVal_tp9,lrt_p=pVal_tp9_lrt)
# 
gene_GE_tp9 <- merge(trans,pvals_tp9_df,by="location")
gene_GE_tp9 <- gene_GE_tp9[order(gene_GE_tp9$wald_p),]
kable(gene_GE_tp9[1:50,]) %>%
  kable_styling()
```

### DEGs specifically for individuals from **TP80**

``` r
ta_tp80 <-  ta_filter[,model$Time == "80"]
dim(ta_tp80)
head(ta_tp80)
model_tp80 <- model[model$Time == "80",]
dim(model_tp80)
head(model_tp80)
tp80 <- DESeqDataSetFromMatrix(countData = ta_tp80,
                              colData = model_tp80,
                              design = ~Lane + Pop + Treatment)

# Wald test using negative binomial
t_tp80_Wald <- DESeq(tp80,minmu = 80,test = "Wald")
# Likelihood ratio test
t_tp80_LTR <- DESeq(tp80,minmu=80,test="LRT",reduced = ~Lane + Pop)

# Treatment at timepoint 9 with Wald
res_t_tp80 <- results(t_tp80_Wald)
summary(res_t_tp80)
pVal_wald_t_tp80 <- res_t_tp80$padj

# Treatment at timepoint 9 with LRT
res_t_tp80_lrt <- results(t_tp80_LTR)
summary(res_t_tp80_lrt)
pVal_lrt_t_tp80 <- res_t_tp80_lrt$padj

# pvals_tp80_df <- data.frame(location=row.names(gc),wald_p=pVal_tp9,lrt_p=pVal_tp9_lrt)
# 
# gene_GE_tp80 <- merge(trans,pvals_tp80_df,by="location")
# gene_GE_tp80 <- gene_GE_tp80[order(gene_GE_tp80$wald_p),]
# kable(gene_GE_tp80[1:50,]) %>%
#   kable_styling()
```
