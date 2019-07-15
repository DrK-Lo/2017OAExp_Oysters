# Markdown files for analyzing oyster 2017 OA exposure experiment 

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
