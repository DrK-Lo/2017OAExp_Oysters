library(parallel)

setwd("/shared_lab/20180226_RNAseq_2017OAExp/DNAm")
#source("scripts/R/DNAmRefCode.R")

sprintf("Reading in Data....")
gene <- readRDS("reference/gene_GeneLoc.RData")
gene <- data.frame(index=1:nrow(gene),gene)
#exon <- readRDS("reference/Exon_GeneLoc.RData")
#exon <- data.frame(index=1:nrow(exon),exon)

meta <- readRDS("processed_samples/03_CytoSummaries/All_CytoSum_summaryTable.RData")
# Create a unique index integer for each cytosine in CytoSummary File
meta$index <- c(1:nrow(meta))
y<-gene
ncores <- 40

# Using chr, position, and strand info for each CpG
# This function matches the corresponding gene

findMatch <- function(chr,pos,strand,motif,tri,index){
  if(any(which(as.character(chr) == y$chr))){
    temp <- y[which(as.character(chr) == y$chr),]
    if(any(which(as.numeric(pos) >= temp$start))){
      temp <- temp[which(as.numeric(pos) >= temp$start),]
      if(any(which(as.character(strand) == temp$strand))){
        temp <- temp[which(as.character(strand) == temp$strand),]
        if(any(which(as.numeric(pos) <= temp$end))){
          temp <- temp[which(as.numeric(pos) <= temp$end),]
          out_temp<-NULL
          for(i in 1:nrow(temp)){
            out_temp <- c(out_temp,c(as.numeric(index),i,temp$index[i]))
          }
          return(out_temp)
          #return(cbind(as.numeric(index),temp$index))
        }else(return(c("NA","NA","NA")))
      }else(return(c("NA","NA","NA")))
    }else(return(c("NA","NA","NA")))
  }else(return(c("NA","NA","NA")))
}

mapply_findMatch <- function(id){
  tmp <- meta[(id[1]+1):id[2],]
  mapply(findMatch,tmp$V1,tmp$V2,tmp$V3,tmp$V6,tmp$V7,tmp$index)
}

id <- floor(
  quantile(0:nrow(meta),
           1-(0:ncores)/ncores
  )
)
idm <- embed(id,2)

sprintf("Starting matching.....")
sprintf(paste0("Working on..",ncores," cores"))

result <-mclapply(nrow(idm):1,
                  function(x) mapply_findMatch(idm[x,]),
                  mc.cores=ncores)
sprintf("Preliminary Index List Saved...")
saveRDS(result,"CytoSummary_Gene_Index_RAWList.RData")

#sprintf("Removing nulls...")
#result <- mclapply(result,
#function(x){x[x!="NA"]},
#function(x){Filter(Negate(is.null),x)},
#mc.cores=ncores)

sprintf("Unlisting and creating matrix...")
index <- matrix(unlist(result),ncol=3,byrow = TRUE)
saveRDS(index,"CytoSummary_Gene_IndexwithNA.RData")

sprintf("Condensing index list (removing NAs)...")
index <- index[index[,1] != "NA",]
saveRDS(index,"CytoSummary_Gene_Index.RData")

sprintf("Subsetting  meta data and gene list with indexes...")
full <- cbind(meta[index[,1],],index[,2],gene[index[,3],])
col<-c("chr_1","pos","cg_strand","motif","tri",
       "cg_index","match_num","gene_index","chr_2","method",
       "feature","start","end","score","gene_strand","phase","gene_ID",
       "gene_id")

sprintf("Saving final data set....")
saveRDS(full,"CytoSummary_Gene_Match.RData")
