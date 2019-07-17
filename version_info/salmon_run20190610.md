# Parameter Info for the analysis of run20190610


## Mapping Info

**Files**

* Base Directory (from project path)
  * `/salmon_files/run20190610_data`
* Raw salmon quantification output (quant.sf files)
  * `/salmon_files/run20190610_data/run20190610`
* Counts (single sample) using `tximport`
  * `/salmon_files/run20190610_data/run20190610_counts`
* `.h5` from `wasabi` (for sleuth)
  * `/salmon_files/run20190610_data/run20190610_h5`

* Mapper: Salmon
* Index : `/salmon_files/oyster_index`

Details for quantification step (output from ` cmd_info.json` from a single sample):
```
{
    "salmon_version": "0.13.1",
    "index": "index/oyster_index",
    "libType": "A",
    "mates1": "/shared_lab/20180226_RNAseq_2017OAExp/raw/17XXX_R1_001.fastq.gz",
    "mates2": "/shared_lab/20180226_RNAseq_2017OAExp/raw/17XXX_R2_001.fastq.gz",
    "threads": "40",
    "validateMappings": [],
    "rangeFactorizationBins": "4",
    "numBootstraps": "1000",
    "gcBias": [],
    "seqBias": [],
    "output": "run20190610/17122",
    "auxDir": "aux_info"
}
```

Example of meta out for sample 17005:
```
{
    "salmon_version": "0.13.1",
    "samp_type": "bootstrap",
    "opt_type": "vb",
    "quant_errors": [],
    "num_libraries": 1,
    "library_types": [
        "ISR"
    ],
    "frag_dist_length": 1001,
    "seq_bias_correct": true,
    "gc_bias_correct": true,
    "num_bias_bins": 4096,
    "mapping_type": "mapping",
    "num_targets": 66451,
    "serialized_eq_classes": false,
    "eq_class_properties": [],
    "length_classes": [
        1082,
	1760,
	2506,
	3859,
	64880
    ],
    "index_seq_hash": "cb0912362f763367efa6bf395643a87cab296cd724b185ce2d9da8eaf9400c9d",
    "index_name_hash": "ce2ed2c35878ac883ca6b767b21893fad5f88f1a7b0fe9bc33ce1678320c8f14",
    "index_seq_hash512": "23a985c5e08222645826d29676463847a3584de7235eb3e5a7fafbe174c4990dba85276b777514a2804428213a0a0b2f3a638c38cba099fb1809f899d50d7cf7",
    "index_name_hash512": "d86153cb9fff9deff4a0015bb79a64d821758bc76ac8c2cac668a7b6abe830b4c24e028a5251ef79a444b2e7faba206d42b3ee2872ae371bdd834e85f85d1ad2",
    "num_bootstraps": 1000,
    "num_processed": 39911842,
    "num_mapped": 16690315,
    "num_dovetail_fragments": 9459768,
    "num_fragments_filtered_vm": 2226335,
    "num_alignments_below_threshold_for_mapped_fragments_vm": 1778105,
    "percent_mapped": 41.817952175697637,
    "call": "quant",
    "start_time": "Mon Jun 10 20:50:07 2019",
    "end_time": "Mon Jun 10 21:19:33 2019"
}
```

### Counts and Matrix Creation (single tximport based script)

Create with script : `scripts/salmon_scripts/R/salmon_countMatrix.R`

* Tximport is used to perform both gene and transcript level counts
  * Primary scaling (for both genes and transcripts) : `countsFromAbundance = "lengthScaledTPM"`

**SessionInfo**
```
R version 3.6.0 (2019-04-26)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS:   /usr/local/programs/R/3.6.0/lib64/R/lib/libRblas.so
LAPACK: /usr/local/programs/R/3.6.0/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] DESeq2_1.24.0               SummarizedExperiment_1.14.0
 [3] DelayedArray_0.10.0         BiocParallel_1.18.0        
 [5] matrixStats_0.54.0          Biobase_2.44.0             
 [7] GenomicRanges_1.36.0        GenomeInfoDb_1.20.0        
 [9] IRanges_2.18.1              S4Vectors_0.22.0           
[11] BiocGenerics_0.30.0         tximport_1.12.1            

loaded via a namespace (and not attached):
 [1] bit64_0.9-7            splines_3.6.0          Formula_1.2-3         
 [4] assertthat_0.2.1       latticeExtra_0.6-28    blob_1.1.1            
 [7] GenomeInfoDbData_1.2.1 pillar_1.4.1           RSQLite_2.1.1         
[10] backports_1.1.4        lattice_0.20-38        glue_1.3.1            
[13] digest_0.6.19          RColorBrewer_1.1-2     XVector_0.24.0        
[16] checkmate_1.9.3        colorspace_1.4-1       htmltools_0.3.6       
[19] Matrix_1.2-17          plyr_1.8.4             XML_3.98-1.20         
[22] pkgconfig_2.0.2        genefilter_1.66.0      zlibbioc_1.30.0       
[25] purrr_0.3.2            xtable_1.8-4           scales_1.0.0          
[28] htmlTable_1.13.1       tibble_2.1.3           annotate_1.62.0       
[31] ggplot2_3.1.1          nnet_7.3-12            lazyeval_0.2.2        
[34] survival_2.44-1.1      magrittr_1.5           crayon_1.3.4          
[37] memoise_1.1.0          foreign_0.8-71         tools_3.6.0           
[40] data.table_1.12.2      stringr_1.4.0          munsell_0.5.0         
[43] locfit_1.5-9.1         cluster_2.0.8          AnnotationDbi_1.46.0  
[46] compiler_3.6.0         rlang_0.3.4            grid_3.6.0            
[49] RCurl_1.95-4.12        rstudioapi_0.10        htmlwidgets_1.3       
[52] bitops_1.0-6           base64enc_0.1-3        gtable_0.3.0          
[55] DBI_1.0.0              R6_2.4.0               gridExtra_2.3         
[58] knitr_1.23             dplyr_0.8.1            bit_1.1-14            
[61] Hmisc_4.2-0            stringi_1.4.3          Rcpp_1.0.1            
[64] geneplotter_1.62.0     rpart_4.1-15           acepack_1.4.1         
[67] tidyselect_0.2.5       xfun_0.7  
```

**Tximport and Bash Script**
```
## Both packages should be installed ahead of time from CRAN
library(tximport)
library(readr)
library(DESeq2)

setwd("/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files")
## NOTE: be sure the both the folder with your salmon transcript quantifications and the file with the transcript summaries are located within this path

### Directory of sample salmon files
DIR <- "run20190610/" # SET THIS to the sample directory with the files you want included in gene count matrix
ls_files <- list.files(DIR)
files_input <- file.path(DIR,ls_files,"/quant.sf",fsep = "")

print("Files being run...")
print(files_input)


### Transcript List
TRANS <- readRDS("NCBI_transcriptome/transcriptome_table.RData") # SET THIS to name of file with your list of transcripts
tr2gene <- data.frame(TXNAME = TRANS$fullID,GENEID=TRANS$location)
tr2gene$TXNAME <- as.character(tr2gene$TXNAME)
tr2gene$GENEID <- as.character(tr2gene$GENEID)

### Sample Metadata (needed for DESeqData Object)
model<-read.csv("metadata/metadata_cvirginica_rna_meta.txt", header=TRUE)
model$Treatment <- as.factor(model$treatment)
model$Time <- as.character(model$timepoint)
model$Time[model$Time == "3"] <- "09"
model$Time[model$Time == "6"] <- "80"
model$Time <- as.factor(model$Time)
model$Pop <- as.factor(model$population)
model$Lane <- as.factor(model$lane)
model$SFV <-interaction(model$Time,model$Treatment) # Creates single factor variable for combination of time and treatment

print("Directory set and transcript list configure...")

### RUN TXIMPORT

print("Starting counting...")

# For gene level quantification
geneAggr <- tximport(files_input,
            type = "salmon",
            countsFromAbundance = "lengthScaledTPM",
            tx2gene = tr2gene)
print("Complete transcript counts...")

saveRDS(tranAggr,"tranMatrixFull.RData",compress = TRUE)
saveRDS(tranAggr$abundance,"tranMatrixAbundance.RData",compress = TRUE)
saveRDS(tranAggr$countsFromAbundance,"tranMatrixCountsFromAbundance.RData",compress = TRUE)
saveRDS(tranAggr$length,"tranMatrixLength.RData",compress = TRUE)

print("Creating DESeqDataSet....")

dds_g <- DESeqDataSetFromTximport(geneAggr,
                                  model,
                                  design = ~ Lane + Pop + Treatment + Time + Treatment:Time)

dds_t <- DESeqDataSetFromTximport(tranAggr,
                                  model,
                                  design = ~ Lane + Pop + Treatment + Time + Treatment:Time)

saveRDS(dds_g,"geneDESeqDataObj.RData")
saveRDS(dd_t,"tranDESeqDataObj.RData")

print("Successfully saved outputs...")
```

## DE Analysis

### DESeq2

### Edgr-Limma

### Sleuth

