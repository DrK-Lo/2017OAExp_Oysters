# Salmon Pipeline

## Overview  
This pipeline takes advantage of a new tool, Salmon, for rapid transcript (and gene) level quantification from RNA-seq data. It creates a standard estimated gene count matrix, which can be used for standard downstream analysis.

### Inputs
1) Raw read data (formated as `.fastq` or `.fastq.gz` or `fq.gz`)
2) Transcriptome File (reference trasncriptome; `.fna` or `.fa` ideally; [NCBI LINK for oyster](https://www.ncbi.nlm.nih.gov/genome/?term=crassostrea+virginica))
3) Gene Annotation File (annotation for each gene; `.gff` or `.gtf` ideally; [NCBI LINK for oyster](https://www.ncbi.nlm.nih.gov/genome/?term=crassostrea+virginica))

### Primary Outputs
1) Tab delimited gene and transcript matrix (matrix, Transcript X # of Samples)
2) `.h5` file for each individual which can be used by sleuth for differential expression.

## Table of Contents 

1. [Brief Description and Literature on Required Tools and Scripts](#one)
2. [Quick Thoughts and Disclaimers before starting](#dis)
3. [Step 1 - Mapping and Transcript Quantification](#two)
    * [Step 1.1: Create and index of potential transcripts](#one.one)
    * [Step 1.2: Perform quasi mapping and quantification](#one.two)
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

**Overview** : Steps takes the raw reads (those after trimming and QC) and partially maps them to a list of reference transcripts (created from the NCBI transcriptome). Next, it probabilistically estimates the abundance of each transcript for each sample.

**Steps** 
 * [Step 1.1: Create and index of potential transcripts](#one.one)
 * [Step 1.2: Perform quasi mapping and quantification](#one.two)

**Required Input**
1) Raw read data (formated as `.fastq` or `.fastq.gz` or `fq.gz`)
2) Transcriptome File (reference trasncriptome; `.gff` or `.gtf` ideally; [NCBI LINK for oyster](https://www.ncbi.nlm.nih.gov/genome/?term=crassostrea+virginica))

**Output**
1) `quant.sf` : Salmons transcript count output format
2) A unique folder for each individual that will also include a number of summary and auxilary files.

#### **Step 1.1 - Create a Index of all transcripts** <a name="one.one"></a>

Create a index of all potential transcripts with the function `salmon index`.

Example Code:
``` 
salmon index -t pathways/and/nameOfTranscriptome.fa.gz -i /pathways/and/nameofFolderForStoringIndex
```

* This will only take a few minutes on a server. 
* Note this only needs to happen once, and you will continue to reference it in the next command for each sample

#### **Step 1.2 - Perform quasi mapping and quantification** <a name="one.two"></a>

Command line for running a Single Sample:
```
> salmon quant -i index -l IU -1 lib_1_1.fq lib_2_1.fq -2 lib_1_2.fq lib_2_2.fq --validateMappings -o out
```

Command line for bash script for running list of samples in a folder:
```
> quant_Salmon.sh /pathway/to/rawfiles /pathway/to/SalmonResultsDirectory /newFolderForResults
```

Script needs three things:
    
    1) A complete pathway to folder with raw files
    2) A pathway to the directory where salmon outputs are stored.
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

### Step 2 - Creating a Gene Count Matrix and `.h5` file <a name="three"></a>

**Overview** : Create a gene count matrix using the R package `tximport` using a custom script. In addition, convert the standard salmon output into an `.h5` formated file so that it can be used by sleuth to perform differential expression analysis.

**Steps** 
 * [Step 2.1: Create a list of transcript and gene levels IDs](#two.one)
 * [Step 2.2: Use `tximport` to created scaled count matrices](#two.two)
 * [Step 2.3: Use `wasabi` to create an `.h5` file](#two.three)

**Step Input**
1) Raw read data (formated as `.fastq` or `.fastq.gz` or `fq.gz`)
2) Transcriptome File (reference trasncriptome; `.fna` or `.fa` ideally; [NCBI LINK for oyster](https://www.ncbi.nlm.nih.gov/genome/?term=crassostrea+virginica))
3) Gene Annotation File (annotation for each gene; `.gff` or `.gtf` ideally; [NCBI LINK for oyster](https://www.ncbi.nlm.nih.gov/genome/?term=crassostrea+virginica))

#### **Step 2.1: Create a list of transcript and gene levels IDs** <a name="two.one"></a>

`tximport` command needs a reference file that can link the transcript ID with the geneID. Here we create a data.frame based on the NCBI reference transcriptome which is formatted appropriately for `tximport`, using a custom script.

[Script Markdown](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/Transcriptomic/RNAseq_Additional/transcriptomeReferenceFile_generation.md)

**2.1 Outputs**
* .RData file with a correctly formated data.frame for linking transcript and gene level IDs


#### **Step 2.2: Use `tximport` to created scaled count matrices** <a name="two.two"></a>

Scaled estimates are created for each sample using the command `tximport` from the packaged `tximport` in R.

**Inputs**
* The `trgene` obj in the gene level estimate is the `.RData` dataframe from the previous step
* `files_input` can be a single or list of `quant.sf` files from a single or group of samples.

**Output**
* R list object that contains five items
    1) `counts` : Counts based on scaling approach 
    2) `abundance` : Abundance 
    3) `length` : Length of each transcript
    4) `countsFromAbundance` : Scaling method 

Example code for **Gene-level** Estimation :
```
geneAggr <- tximport(files_input,
            type = "salmon",
            countsFromAbundance = "lengthScaledTPM",
            tx2gene = tr2gene)
```
* `countsFromAbundance = "lengthScaledTPM"` : Scales coutns using the average transcript length over samples and then the library size (this works for gene level quantification).

Example code for **Transcript-level** Estimation : 
```
tranAggr <- tximport(files_input,
         type = "salmon",
         countsFromAbundance = "dtuScaledTPM",
         txOut = TRUE,
	 tx2gene = tr2gene)
```

* `txOut = TRUE` : indicates you want transcript level adundance estimates
* `countsFromAbundance = "dtuScaledTPM"` : scaled using the median transcript length among isoforms of a gene, and then the library size

[Full R script]()


**NOTE** :  This script also will take your `tximport` list and convert it into a `DESeq obj` and store the output. This can also be created later, but may be convienent since the `tximport` object is large and may not be easily manipulated on a local machine (by comparison the DESeq object is much smaller).

#### **Step 2.3: Use `wasabi` to create an `.h5` file** <a name="two.three"></a>

To perform differential expression using `sleuth` with the outputs from `Salmon`, `wasabi` was used to create a `.h5` object from the `quant.sf` file.

R script
```
library(wasabi)
library(dplyr)

#Using Wasabi to convert Salmon files to ```.h5``` file format for sleuth

dir <- "/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files"
sfdirs <- filepaths("/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files/run20190610",list.files("/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files/run20190610"))
prepare_fish_for_sleuth(sfdirs)
```
[Full R script]()

* **NOTE** : Script also contains additional code for performing sleuth analysis, which require a metadata file that contains sample information.


