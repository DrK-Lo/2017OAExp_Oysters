# Markdown files for analyzing oyster RNA data  

## Folders

* **Generally ordered by position in RNAseq workflow**
* Letters, `A`, indicate differential analyses at the same step.
  
**`1A_STAR_pipeline/`**:  This folder contains scripts for a pipeline to look at differential expression in RNA seq data using a work flow that utilizes:  
- ```Trimmomatic``` implemented using ```dDocent``` to trim adapters and for quality control  
- ```STAR``` to map reads to reference genome  
- ```RSEM``` to determine the the count for each gene for each individual based on aligned reads   
- `outdated/`:  This folder contains previous versions of the bioinformatics pipeline and preliminary analysis. The bioinformatics script has since been broken in to several scripts (listed below) to improve clarity and allow for great pipeline flexibility.

**`1B_Salmon_pipeline/`**:  This folder contains scripts for a pipeline to look at differential expression in RNA seq data using a work flow that utilizes:
- ```Salmon``` to map reads to quantify transcripts using reference transcriptome  
- ```tximport``` pack in R to aggregate transcripts into gene counts.
  
**2A_Multivariate**

**2B_Differential Expression**
- `EdgR and LIMMA` to help with final filtering of gene count matrix, normalize the expression data, and to determine differential expression.
- `DESeq2` to determine differential expression. 


**`3_Phenotype_Analysis/`**: This folder contains the scripts for analyzing the phenotype data.
- `AE_epf_Phenotype.Rmd`: Markdown for the extra-pallial fluid (EPF) pH analysis. 

**`extra/`**: Folder contains scripts ancillary to the gene expression pipeline. Example, a tutorial for set your your `$PATH` variable.
  
## General File structure for STAR Pipeline
  
