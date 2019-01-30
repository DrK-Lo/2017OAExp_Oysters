# Markdown files for analyzing oyster RNA data  

## Folders  
**```outdated\/```**:  This Folder contains previous versions of the bioinformatics pipeline and preliminary analysis. The bioinformatics script has since been broken into several scripts (listed below) to improve clarity and allow for great pipeline flexibility.  
  
## Files  
**```CV17_RNA_trimmingAndAlignment.Rmd```**: File reads in raw RNA seq data and assembled reference genome and performs bot the trimming and quality control steps, as well align remaining reads to the reference. This is primarily accomplished with the alignment pipeline ```dDocent```.  
**```CV17_RNA_refGenome.Rmd```**: Walks through downloading and assembling the oyster reference genome from NCBI.
**```CV17_RNA_countMatrix.Rmd```**: Take the aligned reads generated from the **```CV17_RNA_trimmingAndAlignment.Rmd```** file converts it into a count matrix.  
**```CV17_RNA_countAnalysis.Rmd```**: A brief work up of the count data (currently looks at global patterns of differential expression and individual loci significance).  