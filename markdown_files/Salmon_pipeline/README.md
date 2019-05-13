# Salmon Pipeline

## Overview  
This pipeline takes advantage of a new tool, Salmon, for rapid transcript (and gene) level quantification from RNA-seq data. The leverages the output of Salmon for use in some basic multivariate visualization (RDA, DAPC) and inference based (PERMANOVA) analyses, as well as individual transcript, gene, and isoform level association tests (implemented in both DESeq2 and Sleuth).

### Literature on Required Tools and Scripts

**Mapping and Transcript Quantification**

*Salmon* -  
* [Website](https://salmon.readthedocs.io/en/latest/salmon.html)  
* [Publication]()

**Formatting and Gene Aggregation Tools**  

*Wasabi* - Converts the converts standard salmon outputs into an ```.h5``` file for each individual. These can be used to analyze transcript level data via Sleuth.
* [Github Page](https://github.com/COMBINE-lab/wasabi)

*Tximport* - 
* [Website](http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)

**Gene, Transcript, and Isoform differential expression**

*DESeq2* - 
* [Website](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)


*Sleuth* - 
* [Website](https://pachterlab.github.io/sleuth/about)
* [Pimental et al. 2017 Nature Methods](https://www.nature.com/articles/nmeth.4324)




