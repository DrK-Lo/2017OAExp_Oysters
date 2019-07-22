# General Pipeline for Differential Expression Analysis 

**Overview**: This is intended to be a generalized pipeline, which takes a count matrix (transcripts x # of samples) and a metadata sheet (sample info along with all explainatory and random factors of interest) and performs differential expression to determine which (if any) genes are differentially expressed due to ocean acidification. This is done using two primary workflows: `DESeq2` and `Edgr` with `limma-voom`. Additionally, there is an alternative workflow being explored that uses `sleuth` and the output from `salmon` (IN PROGRESS). 

**Required Inputs**:
* Count matrix (multiple file formats)
    * `.RData` from `tximport` and salmon output
    * tab delimited `.txt` file from `STAR` and `RSEM` output

## Table of Contents 

1. [Brief Description and Literature on Required Tools and Scripts](#dis)
2. [Quick Thoughts and Disclaimers before starting](#qt)
3. [Step 1 - Differential Expression with DESeq2](#one)
4. [Step 2 - Differential Expression with Edgr Limma Voom ](#two)
5. [Step 3 - Examining overlap between workflows](#three)
6. [Step 4 - Target gene differential expression](#four)

### Brief Description and Literature on Required Tools and Scripts <a name="dis"></a>

*DESeq2* - Approach that tests for differential expression using negative binomial generalized linear models. Capable of handling interactions and various types of shrinkage estimators.
* [Website](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
* [Love et al. 2014 Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)

*Edgr*

*limma*

*sleuth* - Aproach that tests for differential expression by taking a model comparison approach, and using likelihood ratio tests to determine the significance of a particular factor relative to a null model (generally one that includes all random or fixed effects not currently being evaluated).
* [Website](https://pachterlab.github.io/sleuth/about)
* [Pimental et al. 2017 Nature Methods](https://www.nature.com/articles/nmeth.4324)

### Quick Thoughts and Disclaimers before starting <a name="qt"></a>
FILLIN

### Step 1 - Differential Expression with DESeq2 <a name="one"></a>

[Script for Salmon Data]()

### Step 2 - Differential Expression with Edgr Limma Voom <a name="two"></a>


### Step 3 - Examining overlap between workflows <a name="three"></a>

### Step 4 -  Step 4 - Target gene differential expression <a name="four"></a>
