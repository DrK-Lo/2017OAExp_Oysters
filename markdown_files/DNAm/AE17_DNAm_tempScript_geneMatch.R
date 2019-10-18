library(parallel)

ncores <- 40

cg <- readRDS("CG_unstranded_summaryTable.RData")
cg$chr<-as.character(cg$chr)
annot <- readRDS("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/reference/gene_exon_combLOC.RData")

findMatch <- function(start,end){
# Identify if CpG is within a gene
    if(any(which(end >= annot$gene_start))){
      annot.red <- annot.red[which(end >= annot$gene_start),]
        if(any(which(start <= annot.red$gene_end))){
          annot.red <- annot.red[which(start <= annot.red$gene_end),]
# Identify if CpGs within genes are in an exon or intronwhi
            if(any(which(end >= annot.red$start))){
              annot.red <- annot.red[which(end >= annot.red$start),]
                if(any(which(start <= annot.red$end))){
                  annot.red <- annot.red[which(start <= annot.red$end),]
                }else(feature <- "Intron")
            }else(feature <- "Intron")
        }else(feature <- "Intergenic")
  }else(feature <- "Intergenic")
  
  data.frame(feature=feature,)
}
mapply_findMatch <- function(range){
  tmp <- cg[(range[1]+1):range[2],]
  mapply(findMatch,tmp$chr,tmp$start,tmp$end)
}

# Divides the total number of cpg loci into even intervals equal to the number 
# of cores you are using.
id <- floor(
  quantile(0:nrow(cg),
           1-(0:ncores)/ncores
  )
)
#Embed transforms the intervales into a ncores X 2 matrix that contains the loci
# interval range for each core (this will be used by mclapply to split the data 
# equally among cores)/
# the rev() just reverses the quantile order to that the ranges 
# are ascending rather than descending in the matrix
idm <- embed(rev(id),2)
idm <- cbind(idm[,2],idm[,1])

sprintf("Starting matching.....")
sprintf(paste0("Working on..",ncores," cores"))

# the paralized lapply function. Here we base the loci ranges to the 
# mapply_findMatch function in order to subset our full data before finding
# matching cases. 
idm[1,]
result <-mclapply(1:nrow(idm),
                  function(x) mapply_findMatch(idm[x,]),
                  mc.cores=ncores)




#### Incomplete code for creating DNAm - gene link file
for(i in 1:5){
  Sys.sleep(0.001)
  print(i)
  temp <- fa[sT_f$chr[i]== fa$chr,]
  temp <- temp[sT_f$end[i]   >= temp$gene_start,]
  temp <- temp[sT_f$start[i] <= temp$gene_end,]
  if(!is.logical(temp)){
    sT_f[i,c(17:24)] <- temp[1,c(12:19)]
  }
}
temp<-c(0,1)
(temp[temp>2]>1)
(nrow(temp[temp>2])>=1)
is.logical(c(0,2))
