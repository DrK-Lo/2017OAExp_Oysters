# Analysis of oyster 2017 OA exposure experiment 

## Folders  

[**`DNAm/`**](https://github.com/epigeneticstoocean/2017OAExp_Oysters/tree/master/markdown_files/DNAm): Folder contains scripts for analyzing the DNA methylation data from the 24 samples on day 9 and 80.

[**`Transcriptomic/`**](https://github.com/epigeneticstoocean/2017OAExp_Oysters/tree/master/markdown_files/Transcriptomic): Folder contains scripts for analyzing the RNAseq data from the 24 samples on day 9 and 80.
- Mapping done with either `STAR` or `Salmon`.
- Multivariate analysis includes PCA and DAPC.
- Differential expression included two pipelines, `DESeq2` and `Edgr` with `limma-voom`.

[**`Phenotype_Analysis/`**](https://github.com/epigeneticstoocean/2017OAExp_Oysters/tree/master/markdown_files/Phenotype_Analysis): This folder contains the scripts for analyzing the phenotype data.
- `AE_epf_Phenotype.Rmd`: R markdown for the extra-pallial fluid (EPF) pH analysis. 

[**`WC/`**](https://github.com/epigeneticstoocean/2017OAExp_Oysters/tree/master/markdown_files/WC): Folder contains the scripts for analyzing the water chemistry data.
- `WC_table.Rmd`: R markdown for creating a table summarizing the water chemistry data.

[**`tutorials/`**](https://github.com/epigeneticstoocean/2017OAExp_Oysters/tree/master/markdown_files/tutorials): Folder contains scripts ancillary to the gene expression pipeline. Example, a tutorial for set your your `$PATH` variable.
