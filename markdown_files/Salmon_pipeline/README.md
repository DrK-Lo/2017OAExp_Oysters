# Salmon Pipeline

## Overview  
This pipeline takes advantage of a new tool, Salmon, for rapid transcript (and gene) level quantification from RNA-seq data. The leverages the output of Salmon for use in some basic multivariate visualization (RDA, DAPC) and inference based (PERMANOVA) analyses, as well as individual transcript, gene, and isoform level association tests (implemented in both DESeq2 and Sleuth).

## Table of Cotents 

1. [Brief Description and Literature on Required Tools and Scripts](#one)
2. [Step 1 - Mapping and Transcript Quantification](#two)
3. [Step 2 - Formating Salmon Outputs,Gene Aggregation, and creating a Transcript to Gene Reference](#three)
4. [Step 3 - Data Visualization and Multivariate Method Analyses](#four)
5. [Step 4 - Differential Expression (transcript, isoform, and gene level)](#five)
6. [Additional Description of Tools](#six)

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

*DESeq2* - Approach that tests for differential expression using negative binomial generalized linear models. Capable of handling interactions and various types of shrinkage estimators.
* [Website](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
* [Love et al. 2014 Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)

*sleuth* - Aproach that tests for differential expression by taking a model comparison approach, and using likelihood ratio tests to determine the significance of a particular factor relative to a null model (generally one that includes all random or fixed effects not currently being evaluated).
* [Website](https://pachterlab.github.io/sleuth/about)
* [Pimental et al. 2017 Nature Methods](https://www.nature.com/articles/nmeth.4324)

### Step 1 - Mapping and Transcript Quantification <a name="two"></a>

**Summary** : Mapping and Quantification is done using the two-phase approach. 
  * Step 1: Use 
  * Step 2: 

**Step 1.1 - Create a Index of all transcripts**


Command line for running a Single Sample:
```
> salmon quant -i index -l IU -1 lib_1_1.fq lib_2_1.fq -2 lib_1_2.fq lib_2_2.fq --validateMappings -o out
```

Bash Script for running multiple samples
```
#!/bin/bash

cd $2;

echo "Saving in directory $2"

for fn in $(ls $1/*R1_001.fastq.gz);
do
samp=$(echo ${fn} | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 1)
echo "Processing Sample ${samp}"
echo "$1/${samp}_R1_001.fastq.gz"

salmon quant -i oyster_index_from_genome \
        -l A \
	-1 $1/${samp}_R1_001.fastq.gz \
        -2 $1/${samp}_R2_001.fastq.gz \
        -p 40 \
        --validateMappings \
        --rangeFactorizationBins 4 \
        --numBootstraps 1000 \
        --gcBias \
        --seqBias \
        -o run20180512/${samp}
done
```

Command line code for running bash script
```

```
### Step 2 - Formating Salmon Outputs,Gene Aggregation, and creating a Transcript to Gene Reference <a name="three"></a>

### Step 3 - Data Visualization and Multivariate Method Analyses <a name="four"></a>

### Step 4 - Differential Expression (transcript, isoform, and gene level) <a name="five"></a>

### Additional Description of Tools <a name="six"></a>





