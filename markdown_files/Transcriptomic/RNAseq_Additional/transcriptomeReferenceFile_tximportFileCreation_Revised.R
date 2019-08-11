library(dplyr)

### Transcript labels based on gnomon only reference ###
transcript <- read.delim("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/transcriptInfo.tab",
                         header = FALSE,skip = 1)

  ### FNA trascriptome from genomic file ###
fa_line <- readLines("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/GCF_002022765.2_C_virginica-3.0_rna_from_genomic.fna")
trans_fg_select <- fa_line[grep(">",fa_line)]
rm(fa_line)
## Split FNA file
#Gene Location
fg_location <- sub(".*gene=(LOC[0-9]+).*", "\\1", trans_fg_select)
# Gene ID
fg_geneID <- sub(".*GeneID:(.+?)].*", "\\1", trans_fg_select)
# Product
fg_product <- sub(".*product=(.+?)].*", "\\1", trans_fg_select)
# Transcript ID
fg_transcriptID <- sub(".*transcript_id=(.+?)].*", "\\1", trans_fg_select)
# Full ID
fg_fullID <- sub(">(.+?)\\s.*", "\\1", trans_fg_select)
# Merge vectors
merged <- data.frame(gene_name = fg_location,
                     product = fg_product,
                     transcript=fg_transcriptID,
                     fullID = fg_fullID)
# Remove duplicate locations
merge_simple <- merged[!duplicated(merged$gene_name), ]

### GTF file (gene annotation) ###
gtf <-  read.delim("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/KM_CV_genome_transcriptome.tab",header=FALSE)
gtf_final <- gtf[as.character(gtf$V2) == "Gnomon" & as.character(gtf$V3) == "transcript",]
rm(gtf)
gtf_final <- gtf_final[,1:11]
colnames(gtf_final) <- c("seqname","source","feature","Start","End","score","strand","frame","transcript_id","gene_id","gene_name")
gtf_final[,1] <- as.character(gtf_final[,1])
gtf_final[,2] <- as.character(gtf_final[,2])
gtf_final[,3] <- as.character(gtf_final[,3])
gtf_final[,9] <- as.character(gtf_final[,9])
gtf_final[,10] <- as.character(gtf_final[,10])
gtf_final[,11] <- as.character(gtf_final[,11])

# Merge GTF data with FNA file (inclusive)
out <- merge(gtf_final,merge_simple,by="gene_name",all.x=TRUE)
error_no_matching <- out[is.na(out$fullID),]
## Save files
saveRDS(out,"/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/STAR_gnomon_gtf.RData")
saveRDS(error_no_matching,"/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/STAR_gnomon_gtf_errorLog.RData")

### Reformat GTF file for the tximport gene file ###
tximportDF <- data.frame(GENEID=out$gene_name,
                         TXNAME=out$transcript_id,
                         gene_id=out$gene_id,
                         chr=out$seqname,start=out$Start,stop=out$End,predict=out$product)
tximportDF$GENEID <- as.character(tximportDF$GENEID)
tximportDF$TXNAME <- as.character(tximportDF$TXNAME)
tximportDF$gene_id <- as.character(tximportDF$gene_id)
tximportDF$chr <- as.character(tximportDF$chr)
tximportDF$predict <- as.character(tximportDF$predict)

saveRDS(tximportDF,"/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/STAR_gnomon_tximportGeneFile.RData")


