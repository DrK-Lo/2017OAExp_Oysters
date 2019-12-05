#a <- matrix(0, 30000, 4)
#install.packages("gplots")
#install.packages("ComplexHeatmap")
library("gplots")
library("superheat")
library(ggplot2)
library(gridExtra)

#a[1:length(a)] <- rnorm(length(a))
#head(a)
#heatmap(a)
setwd("~/Desktop/DNAmethInteractions/")
meta<- readRDS("data/mainEffects_SignificantMatrix.RData")
a <- readRDS("data/MeanMethylationByTreatment_Genic_5xCoveragePerSample.Rdata")


meta$locus <- rownames(a)[meta$index]
head(meta)
str(meta)
meta$ID=as.character(meta$ID)
identical(meta$locus, meta$ID)
sum(meta$locus!=meta$ID)
head(meta)
dim(meta)
head(a) 
dim(a)

sum(duplicated(colnames(meta))) # should be 0
sum(!complete.cases(meta)) # should be 0
sum(duplicated(colnames(a))) # should be 0
sum(!complete.cases(a)) # should be 0
hist(as.matrix(a))

scale_a <- a-rowMeans(a)
head(a)
head(scale_a)
scale_a <- scale_a[meta$index,]
dim(scale_a)

aclust <- hclust(dist(scale_a))
str(aclust)

ng <- 7
groups <- cutree(aclust, k = ng)
(groupnum<- table(groups))

scale_a_info <- data.frame(scale_a, meta, groups)
head(scale_a_info[order(scale_a_info$groups,decreasing = TRUE),])

scale_a_info$Sig <- 0
scale_a_info$Sig[scale_a_info$Time==1 & scale_a_info$Trt==0 & scale_a_info$Interaction==0] <- 1
scale_a_info$Sig[scale_a_info$Time==0 & scale_a_info$Trt==1 & scale_a_info$Interaction==0] <- 2
scale_a_info$Sig[scale_a_info$Time==1 & scale_a_info$Trt==1 & scale_a_info$Interaction==0] <- 3
scale_a_info$Sig[scale_a_info$Interaction==1] <- 4
table(scale_a_info$Sig)
sum(table(scale_a_info$Sig))
dim(scale_a_info)

### Get info
scale_a_info$location1 <- as.numeric(gsub("\\w+_\\w+.\\d_(\\d+)$","\\1",scale_a_info$ID,perl=TRUE))
scale_a_info$chrom <- as.factor(gsub("(\\w+_\\w+.\\d)_\\d+$","\\1",scale_a_info$ID,perl=TRUE))
as.numeric(tapply(scale_a_info$location1, scale_a_info$chrom , max))

sum(is.na(scale_a_info$location1))


for (i in 1:nlevels(scale_a_info$chrom)){
  c <- levels(scale_a_info$chrom)[i]
  wb <- which(scale_a_info$chrom==c)
  new <- scale_a_info$location1[wb]/max(scale_a_info$location1[wb])
  scale_a_info$location2
}


head(scale_a_info$chrom)
M <- list()

for (i in 1:ng){
  print(i)
  png(paste0("InteractionMeanMethGroup",i,".png"), width=8, height=10, units = "in", res=500)
  superheat(as.matrix(scale_a[groups==i,]), 
            row.dendrogram = FALSE,
            heat.lim=c(-0.5,0.5), 
            heat.pal = c("blue", "white","red"),
            heat.pal.values = c(0,0.5,1),
            legend=FALSE,
            title=paste("Group",i,"Num Loci =",groupnum[i]),
            yr = scale_a_info$Sig[groups==i],
            yr.axis.name = "Sig. Code",
            yr.num.ticks = 4,
            yr.point.alpha = 0.5,
            yr.point.size = 0.5
  )
  if (i==1){
    superheat(as.matrix(scale_a[groups==i,]), 
              row.dendrogram = FALSE,
              heat.lim=c(-0.5,0.5), 
              heat.pal = c("blue", "white","red"),
              heat.pal.values = c(0,0.5,1),
              legend=TRUE,
              title=paste("Group",i,"Num Loci =",groupnum[i]),
              yr = scale_a_info$Sig[groups==i],
              yr.axis.name = "Sig. Code",
              yr.num.ticks = 4
  }
  dev.off()
}


superheat(as.matrix(scale_a[groups==i,]), 
          row.dendrogram = FALSE,
          heat.lim=c(-0.5,0.5), 
          heat.pal = c("blue", "white","red"),
          heat.pal.values = c(0,0.5,1),
          legend=TRUE,
          title=paste("Group",i,"Num Loci =",groupnum[i]),
          yr = scale_a_info$location1[groups==i],
          yr.obs.col = as.numeric(scale_a_info$chrom[groups==i]),
          yr.axis.name = "Chrom",
          yr.num.ticks = 4)










# heatmap.2(as.matrix(scale_a[groups==1,]), scale = "none", 
#           col = bluered(length(keytics)-1), breaks=keytics,
#           trace = "none", density.info = "none",
#           dendrogram = "row", labRow = "",
#           Colv=FALSE, key.title="Scaled Methylation",
#           key.xlab="Scaled Methylation", key=TRUE,
#           keysize=1 ,
#           #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
#           key.par=list(mar=c(3.5,0,3,0))
#           # lmat -- added 2 lattice sections (5 and 6) for padding
#           #lmat=rbind(c(5, 4, 2), c(6, 1, 3)), 
#           #lhei=c(2.5, 5), lwid=c(1, 10, 1)
# )
# 
# heatmap.2(as.matrix(scale_a[groups==2,]), scale = "none", 
#           col = bluered(100), 
#           trace = "none", density.info = "none",
#           dendrogram = "row", labRow = "",
#           Colv=FALSE, key.title="Scaled Methylation",
#           key.xlab="Scaled Methylation")
# heatmap.2(as.matrix(scale_a[groups==3,]), scale = "none", 
#           col = bluered(100), 
#           trace = "none", density.info = "none",
#           dendrogram = "row", labRow = "",
#           Colv=FALSE, key.title="Scaled Methylation",
#           key.xlab="Scaled Methylation")
# 
# 
# 
# 
# heatmap.2(as.matrix(a), scale = "none", col = bluered(100), 
#           trace = "both", density.info = "none",
#           dendrogram = "row", labRow = "",
#           Colv=FALSE, key.title="Methylation",
#           key.xlab="Methylation")
# dev.off()
# 
# png("MeanMeth.png", width=10, height=35, units = "in", res=500)
# b2<- heatmap.2(as.matrix(a), scale = "row", col = bluered(100), 
#           trace = "none", density.info = "none",
#           dendrogram = "both", labRow = "",
#           Colv=TRUE, key.title="Scaled\nMethylation",
#           key.xlab="Methylation")
# dev.off()
# 
# ?heatmap
# library(superheat)
