# Salmon Pipeline

## Overview  
This pipeline takes advantage of a new tool, Salmon, for rapid transcript (and gene) level quantification from RNA-seq data. The leverages the output of Salmon for use in some basic multivariate visualization (RDA, DAPC) and inference based (PERMANOVA) analyses, as well as individual transcript, gene, and isoform level association tests (implemented in both DESeq2 and Sleuth).

## Table of Cotents 

1. [Brief Description and Literature on Required Tools and Scripts](#one)
2. [Step 1 - Mapping and Transcript Quantification](#paragraph1)
3. [Step 2 - Formating Salmon Outputs,Gene Aggregation, and creating a Transcript to Gene Reference](#paragraph2)
4. [Step 3 - Data Visualization and Multivariate Method Analyses]()
5. [Step 4 - Differential Expression (transcript, isoform, and gene level)]()
6. [Additional Description of Tools]()

### Brief Description and Literature on Required Tools and Scripts <a name="one"></a>

**Mapping and Transcript Quantification**

*Salmon* -  Fast transcript quantification tool that utilizes a two phase mapping approach, which aligns RNA-seq reads to the a library of target transcripts, called an index, which can be constructed *de novo* or can simply be a published reference transcriptome. 
* [Website](https://salmon.readthedocs.io/en/latest/salmon.html)  
* [Publication](https://www.nature.com/articles/nmeth.4197)

**Formatting and Gene Aggregation Tools**  

*wasabi* - Converts the converts standard salmon outputs into an ```.h5``` file for each individual. These can be used to analyze transcript level data via Sleuth.
* [Github Page](https://github.com/COMBINE-lab/wasabi)

*tximport* - Tool for importing and summarizing transcript-level abundance estimates for transcript and gene level analysis. 
* [Website](http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)
* [Related Paper](https://f1000research.com/articles/4-1521/v1)

**Gene, Transcript, and Isoform differential expression** Brief Description and Literature on Required Tools and Scripts

*DESeq2* - 
* [Website](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

*sleuth* - 
* [Website](https://pachterlab.github.io/sleuth/about)
* [Pimental et al. 2017 Nature Methods](https://www.nature.com/articles/nmeth.4324)

### Step 1 - Mapping and Transcript Quantification

### Step 2 - Formating Salmon Outputs,Gene Aggregation, and creating a Transcript to Gene Reference

### Step 3 - Data Visualization and Multivariate Method Analyses

### Step 4 - Differential Expression (transcript, isoform, and gene level)

### Additional Description of Tools





