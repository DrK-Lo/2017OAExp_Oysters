---
title: "02 Salmon Pipeline"
author: "adowneywall"
date: "5/12/2019"
output: 
  html_document:
    keep_md: true
editor_options: 
  chunk_output_type: console
---



## **Script Description**  

**Brief Overview** : This script performs a series of basic multivariate approaches to explore the whole genome expression profiles of the 24 oyster RNAseq samples. In partiular this includes:
* A permanova (implemented in the package vegan, called adonis)  
* An RDA (implemented in vegan)
* A DAPC with a focus on treatment(implemented in adegenet)

## **Data**  

Steps:  
* Read in dataframe with metadata for each samples including, treatment, time, population, lane of sequencing, and variable that contains each unique combination level between treatment and time.  
* Read in gene (or transcript) abundance matrix generated from the '01_salmon_countMatrix.Rmd' script.  
* Remove uninformative genes from count matrix.  


```r
## Meta Data for the Model
model<-read.csv("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/metadata_cvirginica_rna_meta.txt", header=TRUE)
model$Treatment <- as.factor(model$treatment)
model$Time <- as.character(model$timepoint)
model$Time[model$Time == "3"] <- "09"
model$Time[model$Time == "6"] <- "80"
model$Time <- as.factor(model$Time)
model$Pop <- as.factor(model$population)
model$Lane <- as.factor(model$lane)
model$SFV <-  interaction(model$Time,model$Treatment) # Creates single factor variable for combination of time and treatment

gc <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/salmon_pipeline/run20180512_gene_abundMatrix_.RData")
gc <- gc[-c(1:205),] #removing some no LOC genes (most trna)
```

## PERMANOVA (implements using adonis from vegan package)

* Testing for statistic differential between treatment, time, and treatment:time  


```r
(out <- adonis(t(gc)~Treatment*Time+Pop+Lane,data=model,permutations = 5000))
```

```
## 
## Call:
## adonis(formula = t(gc) ~ Treatment * Time + Pop + Lane, data = model,      permutations = 5000) 
## 
## Permutation: free
## Number of permutations: 5000
## 
## Terms added sequentially (first to last)
## 
##                Df SumsOfSqs  MeanSqs F.Model      R2  Pr(>F)    
## Treatment       1   0.03170 0.031699 1.45160 0.05749 0.02300 *  
## Time            1   0.04832 0.048325 2.21297 0.08764 0.00020 ***
## Pop             2   0.05215 0.026077 1.19415 0.09459 0.09878 .  
## Lane            1   0.02768 0.027681 1.26761 0.05020 0.09738 .  
## Treatment:Time  1   0.02029 0.020291 0.92921 0.03680 0.60968    
## Residuals      17   0.37123 0.021837         0.67328            
## Total          23   0.55138                  1.00000            
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Plotting data with RDA  

### Plot in multivariate space with RDA (treatment and time)

```r
prin_comp<-rda(t(gc), scale=FALSE)
sum_pri <- summary(prin_comp)
pca_scores<-scores(prin_comp)

pca <- prcomp(t(gc))
eigs <- pca$sdev^2

head(sum_pri$species)
```

```
##                        PC1           PC2          PC3           PC4
## LOC111099029 -2.432954e-03  0.0000231773 2.667991e-03  2.121259e-03
## LOC111099030 -5.100276e-01 -0.4343227395 8.016061e-01 -1.232862e+00
## LOC111099031  1.177123e-04  0.0001933143 1.251508e-04  1.356039e-04
## LOC111099032  5.744518e-05  0.0001395553 4.504468e-05  1.122106e-05
## LOC111099033 -5.265188e-02 -0.0214166032 2.029884e-02 -5.044111e-02
## LOC111099034 -1.775297e-03  0.0031607961 4.987734e-03  2.291729e-03
##                        PC5           PC6
## LOC111099029 -1.594811e-03  1.826297e-04
## LOC111099030  1.827979e-01 -4.438335e-01
## LOC111099031  1.059510e-04 -2.218093e-04
## LOC111099032 -3.720566e-06 -3.173615e-05
## LOC111099033 -2.462102e-02 -8.814008e-03
## LOC111099034  2.640098e-03 -5.490798e-04
```

```r
color_comb <- c("deepskyblue2","blue4","firebrick1","darkred") # colors for population 
model$colors <- "" 
model$colors[model$SFV == unique(model$SFV)[1]] <-  color_comb[2]
model$colors[model$SFV == unique(model$SFV)[2]] <-  color_comb[1]
model$colors[model$SFV == unique(model$SFV)[3]] <-  color_comb[4]
model$colors[model$SFV == unique(model$SFV)[4]] <-  color_comb[3]

ordiplot(prin_comp,type="n",
         xlab=paste0("PC1 ",round(eigs[1] / sum(eigs)*100,1),"% Variance Explained)"),
         ylab=paste0("PC2 ",round(eigs[2] / sum(eigs)*100,1),"% Variance Explained)"))
orditorp(prin_comp,display="sites",labels = FALSE,col=model$colors,cex = 2,pch = 16,)
#ordiellipse(prin_comp,model$SFV,conf=0.90,col = color_comb,lwd = 3)
ordispider(prin_comp,model$SFV,col = color_comb,lwd=2.5)
legend(x=100,y=150,legend = c("Day09_Control","Day80_Control","Day09_OA","Day80_OA"),pch = 16,col=color_comb,xpd = .25)
text(x = -538 ,y = -300, paste0("Adonis P_Treatment = ",out$aov.tab$`Pr(>F)`[1],"*"))
text(x = -550 ,y = -360, paste0("Adonis P_Duration = ",out$aov.tab$`Pr(>F)`[2],"*"))
text(x = -540 ,y = -420, paste0("Adonis P_Interaction = ",out$aov.tab$`Pr(>F)`[5]))
```

![](02_salmon_multiVariateVisualization_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

### Plot colored by population  

```r
ordiplot(prin_comp,type="n",
         xlab=paste0("PC1 ",round(eigs[1] / sum(eigs)*100,1),"% Variance Explained)"),
         ylab=paste0("PC2 ",round(eigs[2] / sum(eigs)*100,1),"% Variance Explained)"))
orditorp(prin_comp,display="sites",labels = FALSE,col=model$Pop,cex = 2,pch = 16)
legend(x=100,y=100,legend = c("Ipswich","Rowley 1","Rowley 2"),pch = 16,col=unique(model$Pop),xpd = .25)
```

![](02_salmon_multiVariateVisualization_files/figure-html/unnamed-chunk-4-1.png)<!-- -->
  
### Cumulative Variance Plot  

```r
out <- prcomp(t(gc))
vars <- apply(out$x, 2, var)  
props <- vars / sum(vars)
cumsum(props)
```

```
##       PC1       PC2       PC3       PC4       PC5       PC6       PC7 
## 0.5087059 0.6360367 0.7242154 0.7884668 0.8190354 0.8430982 0.8643221 
##       PC8       PC9      PC10      PC11      PC12      PC13      PC14 
## 0.8829884 0.8971544 0.9107303 0.9227811 0.9343326 0.9450818 0.9542035 
##      PC15      PC16      PC17      PC18      PC19      PC20      PC21 
## 0.9615578 0.9685951 0.9743923 0.9796035 0.9846414 0.9893100 0.9934340 
##      PC22      PC23      PC24 
## 0.9973108 1.0000000 1.0000000
```

```r
plot(cumsum(props)~c(1:24))
```

![](02_salmon_multiVariateVisualization_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

### TWO STEP DAPC: first create discriminant function from TP 9 samples and predict coordinates on df for day 80 samples.  
  
**Creating DF by treatment with first timepoint**  

```r
early_time_counts <- gc[,model$Day == 9]
early_time_meta <- model[model$Day == 9,]

dapc_treatment_10<-dapc(t(early_time_counts),early_time_meta$treatment,n.pca=8,n.da=2)
# PCs = 8
# clusters = 1
early_time_meta$discriminant_treatment_10 <- dapc_treatment_10$ind.coord

ggplot(early_time_meta,aes(discriminant_treatment_10,fill=as.factor(treatment),colour=as.factor(treatment))) + 
  geom_density(alpha=0.1) + xlim(-8,8) + 
  labs(title="Discriminant Function for Treatment on Day 9 (based on 8 PCs)",
       x="Discriminant function 1",
       colour="Treatment",
       fill="Treatment") +
  theme_bw() +
  scale_color_manual(values=c("deepskyblue2","firebrick1")) +
  scale_fill_manual(values=c("deepskyblue2","firebrick1"))
```

![](02_salmon_multiVariateVisualization_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

**Mapping Day 80 samples**  

```r
late_time_counts <- gc[,model$Day == 80]
late_time_meta <- model[model$Day == 80,]

predict_values <- predict.dapc(dapc_treatment_10,t(late_time_counts))
late_time_meta$discriminant_treatment_10 <-predict_values$ind.scores

whole_meta<- rbind(early_time_meta,late_time_meta)

ggplot(whole_meta,aes(discriminant_treatment_10,fill=as.factor(interaction(Day,treatment)),colour=as.factor(interaction(Day,treatment)))) + 
  geom_density(alpha=0.1) + xlim(-8,8) + 
  labs(title="Discriminant Function for Treatment on Day 9 - Mapped with Day 80 Samples",
       x="Discriminant function 1",
       colour="Day.Treatment",
       fill="Day.Treatment") +
  theme_bw() +
  scale_color_manual(values=c("deepskyblue2","blue4","firebrick1","darkred")) +
  scale_fill_manual(values=c("deepskyblue2","blue4","firebrick1","darkred"))
```

![](02_salmon_multiVariateVisualization_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

### DAPC based on combined time*treatment factor  
  
**Creating DF by treatment with first timepoint**  

```r
dapc_SFV_10<-dapc(t(gc),model$SFV,n.pca=8,n.da=3)
# PCs = 10
# clusters = 3
output <- data.frame(Trt=model$Treatment,Time=model$Time,dapc_SFV_10$ind.coord)

ggplot(output,aes(x=LD1,y=LD2,fill=as.factor(interaction(Trt,Time)),colour=as.factor(interaction(Trt,Time)))) + 
  geom_point(aes(size=5)) + #geom_density(alpha=0.1) + #xlim(-28,28) + 
  labs(title="DAPC for Treatment*Time Combination",
       x="Discriminant function 1",
       y="Discriminant function 2",
       colour="Treatment",
       fill="Treatment") +
    theme_bw() +
  scale_color_manual(values=c("deepskyblue2","blue4","firebrick1","darkred")) +
  scale_fill_manual(values=c("deepskyblue2","blue4","firebrick1","darkred"))
```

![](02_salmon_multiVariateVisualization_files/figure-html/unnamed-chunk-8-1.png)<!-- -->
