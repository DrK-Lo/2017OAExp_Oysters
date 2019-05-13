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




```r
model<-read.csv("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/metadata_cvirginica_rna_meta.txt", header=TRUE)
model$Treatment <- as.factor(model$treatment)
model$Time <- as.character(model$timepoint)
model$Time[model$Time == "3"] <- "09"
model$Time[model$Time == "6"] <- "80"
model$Time <- as.factor(model$Time)
model$Pop <- as.factor(model$population)
model$Lane <- as.factor(model$lane)
model$SFV <-  interaction(model$Time,model$Treatment) # Creates single factor variable for combination of time and treatment

gc <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/salmon_RNA/run20180512_gene_countMatrix_.RData")
gc <- gc[-c(1:205),] #removing some no LOC genes (most trna)
```


```r
(out <- adonis(t(gc)~Treatment*Time+Pop+Lane,data=model,permutations = 999))
```

```
## 
## Call:
## adonis(formula = t(gc) ~ Treatment * Time + Pop + Lane, data = model,      permutations = 999) 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##                Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
## Treatment       1   0.03690 0.036898 2.21311 0.08615  0.003 **
## Time            1   0.03884 0.038845 2.32990 0.09070  0.004 **
## Pop             2   0.03303 0.016517 0.99066 0.07713  0.461   
## Lane            1   0.01820 0.018199 1.09155 0.04249  0.254   
## Treatment:Time  1   0.01787 0.017868 1.07175 0.04172  0.281   
## Residuals      17   0.28343 0.016672         0.66180          
## Total          23   0.42827                  1.00000          
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
prin_comp<-rda(t(gc), scale=FALSE)
pca_scores<-scores(prin_comp)
color_comb <- c("deepskyblue2","blue4","firebrick1","darkred") # colors for population 
model$colors <- "" 
model$colors[model$SFV == unique(model$SFV)[1]] <-  color_comb[2]
model$colors[model$SFV == unique(model$SFV)[2]] <-  color_comb[1]
model$colors[model$SFV == unique(model$SFV)[3]] <-  color_comb[4]
model$colors[model$SFV == unique(model$SFV)[4]] <-  color_comb[3]

ordiplot(prin_comp,type="n")
orditorp(prin_comp,display="sites",labels = FALSE,col=model$colors,cex = 2,pch = 16)
#ordiellipse(prin_comp,model$SFV,conf=0.90,col = color_comb,lwd = 3)
ordispider(prin_comp,model$SFV,col = color_comb,lwd=2.5)
legend(x=220,y=-150,legend = c("Day09_Control","Day80_Control","Day09_OA","Day80_OA"),pch = 16,col=color_comb,xpd = .25)
text(x = -538 ,y = -300, paste0("Adonis P_Treatment = ",out$aov.tab$`Pr(>F)`[1],"*"))
text(x = -550 ,y = -360, paste0("Adonis P_Duration = ",out$aov.tab$`Pr(>F)`[2],"*"))
text(x = -540 ,y = -420, paste0("Adonis P_Interaction = ",out$aov.tab$`Pr(>F)`[5]))
```

![](02_salmon_multiVariateVisualization_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

### Cumulative Variance Plot

```r
out <- prcomp(t(gc))
vars <- apply(out$x, 2, var)  
props <- vars / sum(vars)
cumsum(props)
```

```
##       PC1       PC2       PC3       PC4       PC5       PC6       PC7 
## 0.2521577 0.3937149 0.4819127 0.5598985 0.6307584 0.6904549 0.7393132 
##       PC8       PC9      PC10      PC11      PC12      PC13      PC14 
## 0.7858669 0.8141289 0.8412478 0.8626239 0.8814945 0.8994473 0.9152450 
##      PC15      PC16      PC17      PC18      PC19      PC20      PC21 
## 0.9298163 0.9417580 0.9520958 0.9619097 0.9705773 0.9790004 0.9869943 
##      PC22      PC23      PC24 
## 0.9944101 1.0000000 1.0000000
```

```r
plot(cumsum(props)~c(1:24))
```

![](02_salmon_multiVariateVisualization_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

### TWO STEP DAPC: first create discriminant function from TP 9 samples and predict coordinates on df for day 80 samples.  
  
**Creating DF by treatment with first timepoint**  

```r
early_time_counts <- gc[,model$Day == 9]
early_time_meta <- model[model$Day == 9,]

dapc_treatment_10<-dapc(t(early_time_counts),early_time_meta$treatment,n.pca=8,n.da=2)
# PCs = 5
# clusters = 1
early_time_meta$discriminant_treatment_10 <- dapc_treatment_10$ind.coord

ggplot(early_time_meta,aes(discriminant_treatment_10,fill=as.factor(treatment),colour=as.factor(treatment))) + 
  geom_density(alpha=0.1) + xlim(-8,8) + 
  labs(title="Discriminant Function for Treatment on Day 9 (based on 8 PCs)",
       x="Discriminant function 1",
       colour="Treatment",
       fill="Treatment")
```

![](02_salmon_multiVariateVisualization_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

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
       fill="Day.Treatment")
```

![](02_salmon_multiVariateVisualization_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

### DAPC based on combined time*treatment factor  
  
**Creating DF by treatment with first timepoint**  

```r
dapc_SFV_10<-dapc(t(gc),model$SFV,n.pca=8,n.da=3)
# PCs = 10
# clusters = 3
output <- data.frame(Trt=model$Treatment,Time=model$Time,dapc_SFV_10$ind.coord)

ggplot(output,aes(x=LD1,y=LD2,fill=as.factor(interaction(Trt,Time)),colour=as.factor(interaction(Trt,Time)))) + 
  geom_point(aes(size=5)) + #geom_density(alpha=0.1) + #xlim(-28,28) + 
  labs(title="Discriminant Function for Treatment on Day 9",
       x="Discriminant function 1",
       colour="Treatment",
       fill="Treatment")
```

![](02_salmon_multiVariateVisualization_files/figure-html/unnamed-chunk-7-1.png)<!-- -->
