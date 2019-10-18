Evaluation of Expression for Biomineralization Genes
================
adowneywall
10/18/2019

**Libraries**

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
par(mfrow=c(2,2))
geneCounts <- readRDS(paste0(wd,"/results/Transcriptomic/gene_preNormalization_DGEListObj.RData"))
gc_log <- log2(cpm(geneCounts$counts))
tranCounts <- readRDS(paste0(wd,"/results/Transcriptomic/transcript_preNormalization_DGEListObj.RData"))
tc_log <- log2(cpm(tranCounts$counts))

geneDiff <- readRDS(paste0(wd,"/results/Transcriptomic/gene_EBayesObj.RData"))
tranDiff <- readRDS(paste0(wd,"/results/transcript_EBayesObj.RData"))  
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

**Target Gene List - Primarily biomineralization associated genes**

At the moment these are not specifically investigated during the diff.
gene analysis (all gene are examined). I mostly just wanted to confirm
they were present and being expressed in my dataset using a simple
histogram (below)

``` r
### Target List ###
tl <- read.csv(paste0(wd,"/input_files/RNA/references/Target_BiomineralizationGenes.csv"))
# 23 unique biomineralization genes 
# 136 gene locations
tl$Location <- as.character(tl$Location)
```

### **Histogram target vs total genes**

![](05_AE17_RNA_biomineralizationGenes_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

### **Histogram target vs total transcripts**

![](05_AE17_RNA_biomineralizationGenes_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

**Between treatment-time level variation in
genes**

``` r
geneCounts_T1C <-  rowMeans(gc_log[,as.character(meta$SFVrn) == "D09.400"])
geneCounts_T1E <-  rowMeans(gc_log[,as.character(meta$SFVrn) == "D09.2800"])
geneCounts_T2C <-  rowMeans(gc_log[,as.character(meta$SFVrn) == "D80.400"])
geneCounts_T2E <-  rowMeans(gc_log[,as.character(meta$SFVrn) == "D80.2800"])
geneMeanSummary <- cbind(geneCounts_T1C,geneCounts_T1E,
                       geneCounts_T2C,geneCounts_T2E)
geneGrandMean <- abs(rowMeans(geneMeanSummary))

cvCalc <- function(x){sd(x)/mean(x)}
geneCV <- abs(apply(geneMeanSummary,1,cvCalc))

geneMean_bio <- geneGrandMean[!is.na(match(rownames(gc_log),tl$Location))]
geneCV_bio <- geneCV[!is.na(match(rownames(gc_log),tl$Location))]

plot(log2(geneCV)~geneGrandMean,col=adjustcolor("red",alpha.f = 0.2),ylim=c(-10,10),
     xlab="Log2CPM Expression",ylab="log2(CV)",main="Gene expression variation among treatment and time levels",xlim=c(-0.1,14))
points(log2(geneCV_bio)~geneMean_bio,lwd=5,
       col=adjustcolor("blue",alpha.f = 0.5),pch=16)
legend(x=8,y=6,legend=c("All Genes","Biomineralization Genes"),pch=c(1,16),
       col=c(adjustcolor("red",alpha.f = 0.8),adjustcolor("blue",alpha.f = 0.5)),
       bty="n")
```

![](05_AE17_RNA_biomineralizationGenes_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

**Between treatment-time level variation in
transcripts**

``` r
tranCounts_T1C <-  rowMeans(tc_log[,as.character(meta$SFVrn) == "D09.400"])
tranCounts_T1E <-  rowMeans(tc_log[,as.character(meta$SFVrn) == "D09.2800"])
tranCounts_T2C <-  rowMeans(tc_log[,as.character(meta$SFVrn) == "D80.400"])
tranCounts_T2E <-  rowMeans(tc_log[,as.character(meta$SFVrn) == "D80.2800"])
tranMeanSummary <- cbind(tranCounts_T1C,tranCounts_T1E,
                       tranCounts_T2C,tranCounts_T2E)
tranGrandMean <- abs(rowMeans(tranMeanSummary))

cvCalc <- function(x){sd(x)/mean(x)}
tranCV <- abs(apply(tranMeanSummary,1,cvCalc))

tranMean_bio <- tranGrandMean[!is.na(match(tranCounts$genes$GENEID,tl$Location))]
tranCV_bio <- tranCV[!is.na(match(tranCounts$genes$GENEID,tl$Location))]

plot(log2(tranCV)~tranGrandMean,col=adjustcolor("red",alpha.f = 0.2),ylim=c(-10,10),
     xlab="Log2CPM Expression",ylab="log2(CV)",main="Transcripts expression variation among treatment and time levels",xlim=c(-0.1,14))
points(log2(tranCV_bio)~tranMean_bio,lwd=5,
       col=adjustcolor("blue",alpha.f = 0.5),pch=16)
legend(x=8,y=6,legend=c("All Transcripts","Biomineralization Transcripts"),pch=c(1,16),
       col=c(adjustcolor("red",alpha.f = 0.8),adjustcolor("blue",alpha.f = 0.5)),
       bty="n")
```

![](05_AE17_RNA_biomineralizationGenes_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->
