

matchRef <- function(m,ref){
  out <- NULL
  out[1:nrow(m)] <- 0
  for(i in 1:nrow(m)){
    out[i] <-  match(4,c(as.numeric(m[i,1] == ref[,1])+as.numeric(m[i,3] == ref[,7]) +as.numeric(m[i,2] >= ref$start)+as.numeric(m[i,2] <= ref$end)))
  }
  ref2 <- rbind(ref,NA)
  out[is.na(out)] <- nrow(ref2)
  new <- ref2[out,]
  return(new)
}

findMatch <- function(chr,pos,strand,motif,tri,index){
  if(any(which(as.character(chr) == y$chr))){
    temp <- y[which(as.character(chr) == y$chr),]
    if(any(which(as.numeric(pos) >= temp$start))){
      temp <- temp[which(as.numeric(pos) >= temp$start),]
      if(any(which(as.character(strand) == temp$strand))){
        temp <- temp[which(as.character(strand) == temp$strand),]
        if(any(which(as.numeric(pos) <= temp$end))){
          temp <- temp[which(as.numeric(pos) <= temp$end),]
          return(cbind(as.numeric(index),temp))
        }
      }
    }
  }
}

mapply_findMatch <- function(id){
  tmp <- meta[(id[1]+1):id[2],]
  mapply(findMatch,tmp$V1,tmp$V2,tmp$V3,tmp$V6,tmp$V7,tmp$index)
}

parallelRunMatch <- function(meta,y,NCORE){
  ncores <- NCORE
  id <- floor(
    quantile(0:nrow(meta),
             1-(0:ncores)/ncores
    )
  )
  idm <- embed(id,2)

  result <-mclapply(nrow(idm):1,
                    function(x) mapply_findMatch(idm[x,]),
                    mc.cores=ncores)
  result <- lapply(result,function(x){Filter(Negate(is.null),x)})
  return(matrix(unlist(result),byrow = TRUE))
}




