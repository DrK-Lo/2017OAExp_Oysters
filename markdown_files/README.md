# Markdown files for analyzing oyster RNA data  

## Folders  
**```outdated/```**:  This folder contains previous versions of the bioinformatics pipeline and preliminary analysis. The bioinformatics script has since been broken in to several scripts (listed below) to improve clarity and allow for great pipeline flexibility.  
  
**```STAR_pipeline/```**:  This folder contains scripts for a pipeline to look at differential expression in RNA seq data using a work flow that utilizes:  
- ```Trimmomatic``` implemented using ```dDocent``` to trim adapters and for quality control  
- ```STAR``` to map reads to reference genome  
- ```HTSEQ``` to determine the the count for each gene for each individual based on aligned reads   
- ```EdgR and LIMMA``` to help with final filtering of gene count matrix, normalize the expression data, and to determine differential expression.  

**```Salmon_pipeline/```**:  This folder contains scripts for a pipeline to look at differential expression in RNA seq data using a work flow that utilizes:
- ```Salmon``` to map reads to quantify transcripts using reference transcriptome  
- ```tximport``` pack in R to aggregate transcripts into gene counts.
- ```DESeq2``` to determine differential expression.  
- ```WCGNA``` to look at correlated gene networks. 
  
**```HISAT_pipeline/```** (not currently implemented or available):  This folder contains scripts based on Erin Roberts Expression Pipeline and utilizes:  
- ```BBtools```  
- ```HISAT2```  
- ```StringTie```  
- ```DESeq2```   

**```Phenotype_Analysis/```**: This folder contains the scripts for analyzing the phenotype data.
- ```AE_epf_Phenotype.Rmd```: Markdown for the extra-pallial fluid (EPF) pH analysis. 

**```extra/```**: Folder contains scripts ancillary to the gene expression pipeline. Example, a tutorial for set your your ```$PATH``` variable.
  
## General File structure for STAR Pipeline
  
**Note that only the file endings will be consistent between pipelines.**  
  
**```*_createRefGenome.Rmd```**: Walks through downloading and assembling the oyster reference genome from NCBI.  
**```*_trimming.Rmd```**: Script that reads in raw RNA seq data and assembled reference genome and performs both the trimming and quality control steps.  
**```*_ReadMapping.Rmd```**: Script that aligns trimmed reads to a reference genome.  
**```*_Filtering.Rmd```**: Filtering step that removes aligned reads that fail a certain quality threshold (MAPQ level)  
**```*_createcountMatrix.Rmd```**: Take the aligned reads that pass QC and generates a count matrix.  
**```*_countAnalysis.Rmd```**: A brief work up of the count data (currently looks at global patterns of differential expression and individual loci significance).  
