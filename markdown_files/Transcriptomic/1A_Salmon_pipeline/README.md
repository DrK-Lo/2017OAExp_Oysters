# Salmon Pipeline

## Overview  
This pipeline takes advantage of a new tool, Salmon, for rapid transcript (and gene) level quantification from RNA-seq data. It creates a standard estimated gene count matrix, which can be used for standard downstream analysis.

### Inputs
1) Raw read data (formated as `.fastq` or `.fastq.gz` or `fq.gz`)
2) Transcriptome File (reference trasncriptome; `.gff` or `.gtf` ideally; [NCBI LINK for oyster](https://www.ncbi.nlm.nih.gov/genome/?term=crassostrea+virginica))

### Primary Outputs
1) Tab delimited gene and transcript files (matrix, Transcript X # of Samples)
2) `.h5` file for each individual which can be used by sleuth for differential expression.

## Table of Contents 

1. [Brief Description and Literature on Required Tools and Scripts](#one)
2. [Quick Thoughts and Disclaimers before starting](#dis)
3. [Step 1 - Mapping and Transcript Quantification](#two)
4. [Step 2 - Formating Salmon Outputs](#three)
5. [Step 3 - Gene Aggregation](#four)
6. [Step 4 - Creating a Transcript to Gene Reference](#five)


### Brief Description and Literature on Required Tools and Scripts <a name="one"></a>

**Mapping and Transcript Quantification**

*Salmon* - Fast transcript quantification tool that utilizes a two phase mapping approach, which aligns RNA-seq reads to the a library of target transcripts, called an index, which can be constructed *de novo* or can simply be a published reference transcriptome. 
* [Website](https://salmon.readthedocs.io/en/latest/salmon.html)  
* [Publication](https://www.nature.com/articles/nmeth.4197)

**Formatting and Gene Aggregation Tools**  

*wasabi* - Converts the standard salmon outputs into an ```.h5``` file for each individual. These can be used to analyze transcript level data via Sleuth.
* [Github Page](https://github.com/COMBINE-lab/wasabi)

*tximport* - Tool for importing and summarizing transcript-level abundance estimates for transcript and gene level analysis. 
* [Website](http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)
* [Related Paper](https://f1000research.com/articles/4-1521/v1)

### Quick Thoughts and Disclaimers before starting <a name="dis"></a>

* Scripts are presented as if primary commands (i.e. `salmon quant`) have been placed in your `$PATH` variable.

* Alternatively, salmon specific commands may be installed via conda [see LINK for details](https://combine-lab.github.io/salmon/getting_started/). Be sure to activate salmon prior to beginning pipeline by entering `conda activate salmon` in the command line.


### Step 1 - Mapping and Transcript Quantification <a name="two"></a>

**Summary** : Mapping and Quantification is done using the two-phase approach. 
 * [Step 1.1: Create and index of potential transcripts](#"one.one")
 * [Step 1.2: Performing quasi mapping and quantification](#"one.two")

#### **Step 1.1 - Create a Index of all transcripts** <a name="one.one"></a>

This is with a simple one line code and will only take a few minutes on a server. Note this only needs to happen once, and you will continue to reference it in the next command for each sample:
``` 
salmon index -t pathways/and/nameOfTranscriptome.fa.gz -i /pathways/and/nameofFolderForStoringIndex
```

#### **Step 1.2 - Performing quasi mapping and quantification** <a name="one.two"></a>

Command line for running a Single Sample:
```
> salmon quant -i index -l IU -1 lib_1_1.fq lib_2_1.fq -2 lib_1_2.fq lib_2_2.fq --validateMappings -o out
```

Command line for bash script for running list of samples in a folder:
```
> quant_Salmon.sh /pathway/to/rawfiles /pathway/to/SalmonResultsDirectory /newFolderForResults
```
Script needs three things
	1) A complete pathway to folder with raw files
	2) A pathway to the directory where salmon outputs are stored
	3) A **new** folder name that will be created in the directory to store outputs that are generated

Bash Script of `quant_Salmon.sh` script
```
#!/bin/bash

cd $2;
mkdir $3;

echo "Creating folder $3 for output"
echo "Saving in directory $2"
echo "Retrieving fastq files from $1 directory"

for fn in $(ls $1/*R1_001.fastq.gz);
do
samp=$(echo ${fn} | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 1)
echo "Processing Sample ${samp}"
echo "$1/${samp}_R1_001.fastq.gz"

salmon quant -i index/oyster_index \
        -l A \
	-1 $1/${samp}_R1_001.fastq.gz \
        -2 $1/${samp}_R2_001.fastq.gz \
        -p 40 \
        --validateMappings \
        --rangeFactorizationBins 4 \
        --numBootstraps 1000 \
        --gcBias \
        --seqBias \
        -o $3/${samp}
done
```

### Step 2 - Creating a Gene Count Matrix <a name="three"></a>

Creating a gene count matrix was done using the R package `tximport`. This was performed using a custom script, which was developed to run on a remote cluster and execute the `tximport` on a list of samples from a directory and assemble them into a single count matrix. The specified directory should have a list of quant.sf formated files, which are required for the tximport function.

### Step 3 - Gene Aggregation <a name="four"></a>

### Step 4 - Transcript to Gene Reference <a name="five"></a>

### Additional Description of Tools <a name="six"></a>





