# General Pipeline for Differential Expression Analysis using three primary methods (Edgr-Limma-Voom, DESeq2, and Sleuth)



**Gene, Transcript, and Isoform differential expression** Brief Description and Literature on Required Tools and Scripts

*DESeq2* - Approach that tests for differential expression using negative binomial generalized linear models. Capable of handling interactions and various types of shrinkage estimators.
* [Website](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
* [Love et al. 2014 Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)

*sleuth* - Aproach that tests for differential expression by taking a model comparison approach, and using likelihood ratio tests to determine the significance of a particular factor relative to a null model (generally one that includes all random or fixed effects not currently being evaluated).
* [Website](https://pachterlab.github.io/sleuth/about)
* [Pimental et al. 2017 Nature Methods](https://www.nature.com/articles/nmeth.4324)
