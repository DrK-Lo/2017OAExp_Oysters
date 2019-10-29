
setwd("/shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_files/run_20190726_RSEMGnomonIndex/m3")
lf <- list.files(pattern="exonCounts")
initial<- read.delim(lf[1],header = FALSE)
summary_stats <- grep(pattern="__",initial$V1)
labels<-as.character(initial[,1])
sample_id <- sub("(.+?)_.*","\\1",lf)

# Empty matrices
countMatrix<-matrix(ncol=length(lf),nrow=nrow(initial)-length(summary_stats))
summaryMatrix<-matrix(ncol=length(lf),nrow=length(summary_stats))
row.names(countMatrix) <- labels[-summary_stats]
row.names(summaryMatrix) <- labels[summary_stats]
colnames(countMatrix) <- sample_id
colnames(summaryMatrix) <- sample_id

for( i in 1:length(lf)){
  temp <- read.delim(lf[i],header = FALSE)
  countMatrix[,i] <- temp[-summary_stats,2]
  summaryMatrix[,i] <- temp[summary_stats,2]
}

saveRDS(summaryMatrix,"exonSummaryStatsMatrix.RData")
saveRDS(countMatrix,"exonCountMatrix.RData")