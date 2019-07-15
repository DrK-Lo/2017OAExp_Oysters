Gene Count Matrix Creation
================
adowneywall
6/8/2019

### User set directory information

``` r
# Pathway to folder with samples to be included in gene count matrix
PATH_IN <- "/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/STAR_pipeline/rawCounts"
#Name of folder with raw count files (from STAR)
FOLD <- "rerun3"
# Patway where created .rds files are stored
PATH_OUT <- "/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/STAR_pipeline/rawCounts"
# GCF file used to creat gene counts in STAR
GCF <- read.delim("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/GCF.gtf")
head(GCF)
```

    ##   NC_035780.1 Gnomon transcript X13578 X14594 . X. ..1
    ## 1 NC_035780.1 Gnomon       exon  13578  13603 .  +   .
    ## 2 NC_035780.1 Gnomon       exon  14237  14290 .  +   .
    ## 3 NC_035780.1 Gnomon       exon  14557  14594 .  +   .
    ## 4 NC_035780.1 Gnomon transcript  28961  33324 .  +   .
    ## 5 NC_035780.1 Gnomon       exon  28961  29073 .  +   .
    ## 6 NC_035780.1 Gnomon       exon  30524  31557 .  +   .
    ##   transcript_id.rna0..gene_id.gene0..gene_name.LOC111116054.
    ## 1 transcript_id rna0; gene_id gene0; gene_name LOC111116054;
    ## 2 transcript_id rna0; gene_id gene0; gene_name LOC111116054;
    ## 3 transcript_id rna0; gene_id gene0; gene_name LOC111116054;
    ## 4 transcript_id rna1; gene_id gene1; gene_name LOC111126949;
    ## 5 transcript_id rna1; gene_id gene1; gene_name LOC111126949;
    ## 6 transcript_id rna1; gene_id gene1; gene_name LOC111126949;

### Code for looping through each file to create a single gene count matrix for each time of count

``` r
files <- list.files(paste0(PATH_IN,"/",FOLD))
LEN <- length(files)
NAME <- c("Gene","unstranded","Stranded","Reverse")
for(j in 2:4){
  for(i in 1:LEN){
    if(i == 1){
      TEMP <- read.delim(paste0(PATH_IN,"/",FOLD,"/",files[i]),
                       header = FALSE,skip = 4)
      names(TEMP) <- NAME
      MAT <- matrix(nrow=length(TEMP$Gene),ncol = LEN)
      MAT[,i] <- TEMP[,j]
    }else{
        TEMP <- read.delim(paste0(PATH_IN,"/",FOLD,"/",files[i]),
                       header = FALSE,skip = 4)
        MAT[,2] <- TEMP[,j]
    }
  }
  saveRDS(MAT,paste0(PATH_OUT,"/",FOLD,"_GeneCountMatrix_",NAME[j],".rds"))
}
```

**WARNING**: these files need to be the readCount tab delimited files
produced by STAR
