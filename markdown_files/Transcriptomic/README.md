# RNAseq data analysis  

## Folders

**Generally ordered by position in RNAseq workflow** Letters, `A`, indicate differential analyses at the same step.

**`1A_Salmon_pipeline/`**:  This folder contains scripts for a pipeline to look at differential expression in RNA seq data using a work flow that utilizes:
- ```Salmon``` to map reads to quantify transcripts using reference transcriptome  
- ```tximport``` pack in R to aggregate transcripts into gene counts.

**`1B_STAR_pipeline/`**:  This folder contains scripts for a pipeline to look at differential expression in RNA seq data using a work flow that utilizes:  
- ```Trimmomatic``` implemented using ```dDocent``` to trim adapters and for quality control  
- ```STAR``` to map reads to reference genome  
- ```RSEM``` to determine the the count for each gene for each individual based on aligned reads   
- `outdated/`:  This folder contains previous versions of the bioinformatics pipeline and preliminary analysis. The bioinformatics script has since been broken in to several scripts (listed below) to improve clarity and allow for great pipeline flexibility.
  
**`2A_Multivariate_Visualization/`**:

**`2B_DiffExp/`**
- `EdgR and LIMMA` to help with final filtering of gene count matrix, normalize the expression data, and to determine differential expression.
- `DESeq2` to determine differential expression. 

**`Diagnostics/`**: Folder contains scripts for aiding in some basic diagnostics and comparisons in the gene expression pipeline.

**`RNAseq_Additional/`**: Some additional helpful scripts.
  
