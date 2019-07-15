---
title: "tximport_example"
author: "adowneywall"
date: "5/12/2019"
output: 
  html_document: 
    keep_md: true
editor_options: 
  chunk_output_type: console
---



Basic Options

```r
dir <- system.file("extdata", package = "tximportData")
list.files(dir)
```

```
##  [1] "cufflinks"               "derivedTxome"           
##  [3] "kallisto"                "kallisto_boot"          
##  [5] "rsem"                    "sailfish"               
##  [7] "salmon"                  "salmon_gibbs"           
##  [9] "samples_extended.txt"    "samples.txt"            
## [11] "tx2gene.csv"             "tx2gene.gencode.v27.csv"
```

```r
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samples
```

```
##   pop center                assay    sample experiment       run
## 1 TSI  UNIGE NA20503.1.M_111124_5 ERS185497  ERX163094 ERR188297
## 2 TSI  UNIGE NA20504.1.M_111124_7 ERS185242  ERX162972 ERR188088
## 3 TSI  UNIGE NA20505.1.M_111124_6 ERS185048  ERX163009 ERR188329
## 4 TSI  UNIGE NA20507.1.M_111124_7 ERS185412  ERX163158 ERR188288
## 5 TSI  UNIGE NA20508.1.M_111124_2 ERS185362  ERX163159 ERR188021
## 6 TSI  UNIGE NA20514.1.M_111124_4 ERS185217  ERX163062 ERR188356
```

```r
files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample", 1:6)
all(file.exists(files))
```

```
## [1] TRUE
```

```r
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
```

```
## Parsed with column specification:
## cols(
##   TXNAME = col_character(),
##   GENEID = col_character()
## )
```

```r
head(tx2gene)
```

```
## # A tibble: 6 x 2
##   TXNAME            GENEID           
##   <chr>             <chr>            
## 1 ENST00000456328.2 ENSG00000223972.5
## 2 ENST00000450305.2 ENSG00000223972.5
## 3 ENST00000473358.1 ENSG00000243485.5
## 4 ENST00000469289.1 ENSG00000243485.5
## 5 ENST00000607096.1 ENSG00000284332.1
## 6 ENST00000606857.1 ENSG00000268020.3
```

```r
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
```

```
## reading in files with read_tsv
```

```
## 1 2 3 4 5 6 
## summarizing abundance
## summarizing counts
## summarizing length
```

```r
names(txi)
```

```
## [1] "abundance"           "counts"              "length"             
## [4] "countsFromAbundance"
```
  
Rerun with gibs sampler options from salmon  

```r
files <- file.path(dir, "salmon_gibbs", samples$run, "quant.sf")
all(file.exists(files))
```

```
## [1] TRUE
```

```r
names(files) <- paste0("sample", 1:6)
txi.inf.rep <- tximport(files, type = "salmon", txOut = TRUE)
```

```
## reading in files with read_tsv
```

```
## 1 2 3 4 5 6
```

```r
names(txi.inf.rep)
```

```
## [1] "abundance"           "counts"              "infReps"            
## [4] "length"              "countsFromAbundance"
```
