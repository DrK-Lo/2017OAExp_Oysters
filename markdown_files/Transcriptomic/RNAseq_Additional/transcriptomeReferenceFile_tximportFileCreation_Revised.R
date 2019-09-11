library(dplyr)

### Transcript labels based on gnomon only reference ###
transcript <- read.delim("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/transcriptInfo.tab",
                         header = FALSE,skip = 1)

### FNA trascriptome from genomic file ###
fa_line <- readLines("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/GCF_002022765.2_C_virginica-3.0_rna_from_genomic.fna")
#fa_line <- readLines("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/C_virginica-3.0_rna.fa")

trans_fg_select <- fa_line[grep(">",fa_line)]
rm(fa_line)
trans_fg_select[50]
## Split FNA file
#Gene Location
fg_location <- gsub(".*gene=(LOC[0-9]+).*", "\\1", trans_fg_select)
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

## Read in GO Slim terms ##
goT <- read.delim("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/XP_sequences_Cvirginica_GCF_002022765.2_GO_Final.tab")



# Remove duplicate locations
merge_simple <- merged[!duplicated(merged$gene_name), ]
merge_order <- merged[order(merged$gene_name),]
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
#saveRDS(out,"/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/STAR_gnomon_gtf.RData")
#saveRDS(error_no_matching,"/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/STAR_gnomon_gtf_errorLog.RData")

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

#saveRDS(tximportDF,"/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/STAR_gnomon_tximportGeneFile.RData")

### Subset genome GFF file and adding GO terms #####
fa_gff <- read.delim("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/KM_002022765.2_CV_3.0_geome.gff",sep = "\t",skip=8,header=FALSE)
## CDS with GOTerms ##
## Creating GO term protein - gene LOC file ##
fa_gff_CDS <- fa_gff[as.character(fa_gff$V3)=="CDS",]
fa_gff_CDS <- fa_gff_CDS[as.character(fa_gff_CDS$V2)=="RefSeq"]
CDS_attr <- sub("(.*protein_id=.{14}?).*","\\1",as.character(fa_gff_CDS$V9),perl=TRUE)
#Subsetting different attributes
gff_ID <- sub("ID=(.+?);Parent=(.+?);.*gene=(.+?);.*protein_id=(.{14})","\\1",CDS_attr,perl=TRUE)
gff_parent <- sub("ID=(.+?);Parent=(.+?);.*gene=(.+?);.*protein_id=(.{14})","\\2",CDS_attr,perl=TRUE)
gff_gene <- sub("ID=(.+?);Parent=(.+?);.*gene=(.+?);.*protein_id=(.{14})","\\3",CDS_attr,perl=TRUE)
gff_Prot <- sub("ID=(.+?);Parent=(.+?);.*gene=(.+?);.*protein_id=(.{14})","\\4",CDS_attr,perl=TRUE)

CDS_full <- data.frame(chr=as.character(fa_gff_CDS$V1),method=as.character(fa_gff_CDS$V2),
                       feature=as.character(fa_gff_CDS$V3),start=as.character(fa_gff_CDS$V4),
                       end=as.character(fa_gff_CDS$V5),score=as.character(fa_gff_CDS$V6),
                       strand=as.character(fa_gff_CDS$V7),phase=as.character(fa_gff_CDS$V8),
                       CDS_ID=gff_ID,Parent=gff_parent,
                       gene_id=gff_gene,protein_id=gff_Prot)
go_file <- read.delim("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/XP_sequences_Cvirginica_GCF_002022765.2_GO_Final.tab",colClasses = c(rep("character",times=6)))
CDS_wGO <- merge(CDS_full,go_file,by.x = "protein_id",by.y = "SeqName",all.x = TRUE)
saveRDS(CDS_wGO,"/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/CDS_wGoTerms_GeneLOC.RData")

## For Exons ##
fa_gff_exon <- fa_gff[as.character(fa_gff$V3)=="exon",]
fa_gff_exon <- fa_gff_exon[as.character(fa_gff_exon$V2)!="RefSeq",]
exon_attr <- sub("(.*transcript_id=.{14}?).*","\\1",as.character(fa_gff_exon$V9),perl=TRUE)
#Subsetting different attributes
gff_ID <- sub("ID=(.+?);Parent=(.+?);.*gene=(.+?);.*transcript_id=(.{14})","\\1",exon_attr,perl=TRUE)
gff_parent <- sub("ID=(.+?);Parent=(.+?);.*gene=(.+?);.*transcript_id=(.{14})","\\2",exon_attr,perl=TRUE)
gff_gene <- sub("ID=(.+?);Parent=(.+?);.*gene=(.+?);.*transcript_id=(.{14})","\\3",exon_attr,perl=TRUE)
gff_Prot <- sub("ID=(.+?);Parent=(.+?);.*gene=(.+?);.*transcript_id=(.{14})","\\4",exon_attr,perl=TRUE)

exon_full <- data.frame(cbind(chr=as.character(fa_gff_exon$V1),
                        method=as.character(fa_gff_exon$V2),
                       feature=as.character(fa_gff_exon$V3),
                       start=as.character(fa_gff_exon$V4),
                       end=as.character(fa_gff_exon$V5),
                       score=as.character(fa_gff_exon$V6),
                       strand=as.character(fa_gff_exon$V7),
                       phase=as.character(fa_gff_exon$V8),
                       exon_ID=as.character(gff_ID),
                       Parent=as.character(gff_parent),
                       gene_id=as.character(gff_gene),
                       transcript_id=as.character(gff_Prot)))

exon_full$chr <- as.character(exon_full$chr)
exon_full$start <- as.integer(as.character(exon_full$start))
exon_full$end <- as.integer(as.character(exon_full$end))
exon_full$strand <- as.character(exon_full$strand)
exon_full$exon_ID <- as.character(exon_full$exon_ID)
exon_full$Parent <- as.character(exon_full$Parent)
exon_full$gene_id <- as.character(exon_full$gene_id)
exon_full$transcript_id <- as.character(exon_full$transcript_id)
saveRDS(exon_full,"/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/Exon_GeneLoc.RData")

## Gene ##
fa_gff_gene <- fa_gff[as.character(fa_gff$V3)=="gene",]
fa_gff_gene <- fa_gff_gene[as.character(fa_gff_gene$V2)!="RefSeq",]
gene_attr <- sub("(.*gene_biotype=.{*}?).*","\\1",as.character(fa_gff_gene$V9),perl=TRUE)
#Subsetting different attributes
gff_ID <- sub("ID=(.+?);.*gene=(.+?);.*","\\1",gene_attr,perl=TRUE)
#gff_parent <- sub("ID=(.+?);Parent=(.+?);.*gene=(.+?);.*gene_biotype=(.{14})","\\2",gene_attr,perl=TRUE)
gff_gene <- sub("ID=(.+?);.*gene=(.+?);.*","\\2",gene_attr,perl=TRUE)
#gff_Prot <- sub("ID=(.+?);.*gene=(.+?);.*gene_biotype=(.{*})","\\3",gene_attr,perl=TRUE)

gene_full <- data.frame(chr=as.character(fa_gff_gene$V1),method=as.character(fa_gff_gene$V2),
                        feature=as.character(fa_gff_gene$V3),start=as.character(fa_gff_gene$V4),
                        end=as.character(fa_gff_gene$V5),score=as.character(fa_gff_gene$V6),
                        strand=as.character(fa_gff_gene$V7),phase=as.character(fa_gff_gene$V8),
                        gene_ID=gff_ID,gene_id=gff_gene)

gene_full$chr <- as.character(gene_full$chr)
gene_full$start <- as.integer(as.character(gene_full$start))
gene_full$end <- as.integer(as.character(gene_full$end))
gene_full$strand <- as.character(gene_full$strand)
gene_full$gene_ID <- as.character(gene_full$gene_ID)
gene_full$gene_id <- as.character(gene_full$gene_id)
saveRDS(gene_full,"/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/gene_GeneLoc.RData")
