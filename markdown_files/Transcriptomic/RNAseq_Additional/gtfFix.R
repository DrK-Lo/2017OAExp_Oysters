
#### SET BEFORE RUNNING ####

### SET WORKING DIRECTORY ###
setwd("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references")
### SET NAME OF FILE ###
fileName <- "KM_CV_genome.gtf"
fn_simple <- strsplit(fileName,".gtf")

### Outputs ###

# This will create a revised versions of the gtf file that include random name inserts for missing gene_ids or gene_names
# Example : gene_id 'geneInsert1'
## Files
# _edit.gtf : full .gtf file with filled in missing ids
# _edit_transcript.gtf : just transcripts
# _edit_exon.gtf : just exons
# _edit_Gnomon.gtf : Just those identified from Gnomon
# _edit_noRefSeq.gtf : With transcripts identified with RefSeq removed.


#### DO NOT MODIFY BELOW THIS ####
## Read in Data ##
line <- read.delim(fileName,header = FALSE)
line2 <- readLines(fileName)

## Create matrix for storing heirarchical gene features names ##
new_mat <- matrix(ncol=8,nrow=length(line$V9))
## split the names colume into a single column for each feature
for(i in 1:length(line$V9)){
  new_mat[i,] <- unlist(strsplit(as.character(line$V9[i]),";| "))[1:8]
}

## Identify and Relabel missing Data ##
# Identify rows with missing gene_id
gene_id_index <- which(is.na(new_mat[,5]))
gene_id <- NULL
# Create automated name for missing gene_id and store in vector
for(i in 1:length(gene_id_index)){
  gene_id[i] <- paste0("geneInsert",i)
}
new_mat[gene_id_index,4] <- "gene_id"
new_mat[gene_id_index,5] <- gene_id
# Create automated name for missing gene_name and store in vector
gene_id_index2 <- which(is.na(new_mat[,8]))
gene_name <- NULL
for(i in 1:length(gene_id_index2)){
  gene_name[i] <- paste0("tempLOC",i)
}
new_mat[gene_id_index2,7] <- "gene_name"
new_mat[gene_id_index2,8] <- gene_name

# Create vector for index values for rows missing either gene_id or gene_name
comb_gene_index <- union(gene_id_index,gene_id_index2)

## Run through all ids and reformat ##
# Transcript_ids reformat
transcript_id <- NULL
for(i in 1:nrow(new_mat)){
  transcript_id[i] <- paste0("transcript_id \"",new_mat[i,2],"_",new_mat[i,5],"\"")
}
# gene_ids reformat
gene_id_rv <- NULL
for(i in 1:nrow(new_mat)){
  gene_id_rv[i] <- paste0(" gene_id \"",new_mat[i,5],"\"")
}
# gene_name reformat
gene_name_rv <- NULL
for(i in 1:nrow(new_mat)){
  gene_name_rv[i] <- paste0(" gene_name \"",new_mat[i,8],"\"")
}

# Combine transcript_id,gene_id, and gene_name into single string
atr_string <- paste0(transcript_id,";",gene_id_rv,";",gene_name_rv,";")

# Split lines from original .gtf file, keeping only the columns prior to the labels (will me concentaneted with new string)
initial <- unlist(strsplit(line2,split = "transcript_id"))[seq(1,length(line2)*2,by=2)]
# Combine initial info from original table with new labels
combine <- paste0(initial,atr_string)

# Save whole revised .gtf
fileConn<-file(paste0(fn_simple,"_edit.gtf"))
writeLines(combine,fileConn)
close(fileConn)

# Save transcripts
combine_transcript <- combine[line$V3 == 'transcript']
fileConn<-file(paste0(fn_simple,"_edit_transcript.gtf"))
writeLines(combine_transcript,fileConn)
close(fileConn)

# Save exons
combine_exon <- combine[line$V3 == 'exon']
fileConn<-file(paste0(fn_simple,"_edit_exon.gtf"))
writeLines(combine_exon,fileConn)
close(fileConn)

# Save only those identified by Gnomon
combine_Gnomon <- combine[line$V2 == 'Gnomon']
fileConn<-file(paste0(fn_simple,"_edit_Gnomon.gtf"))
writeLines(combine_Gnomon,fileConn)
close(fileConn)

# Save only not identified by RefSeq
combine_noRefSeq <- combine[line$V2 != 'RefSeq']
fileConn<-file(paste0(fn_simple,"_edit_noRefSeq.gtf"))
writeLines(combine_noRefSeq,fileConn)
close(fileConn)

# Save everything except Genomom
combine_noGnomonSeq <- combine[line$V2 != 'Gnomon']
fileConn<-file(paste0(fn_simple,"_edit_noGnomon.gtf"))
writeLines(combine_noGnomonSeq,fileConn)
close(fileConn)
length(line$V3[line$V3 == "transcript"])
length(unique(new_mat[,2]))
length(new_mat[,2])
