# Pipeline for salmon run : `run20190610`

### Create new index based on latest NCBI transcriptome annotations

Created a new index from the most recent oyster transcriptome (NOT from genome): `GCF_002022765.2_C_virginica-3.0_rna.fna.gz`

Code for creating index:

```
salmon index -t NCBI*/*rna.fna* -i index/oyster_index 
```

* NOTE 1: This code was run from the working directory `/shared*/2018*/RNA/salmon*` on comp5 computer cluster.
* NOTE 2: You made need "activate" salmon within conda using 'conda activate salmon` prior to running script.

### Transcript Quantification

Perform Transcript quantification using `salmon_quant.sh` script (complete code below).

Command line script:
```
salmon_quant /shared*/2018*/raw /shared*/2018*/RNA/salmon* run20190610
```

`/shared*/2018*/raw` : location of raw fastq read files (these are not ones that have been trimmed for use in STAR)

`/shared*/2018*/RNA/salmon*` : path to salmon directory

`run20190610` : new directory name for output (located in based directory)

* NOTE 1: salmon_quant is the name of the alias I created for running the `salmon_quant.sh` script 

Bash Code Script for salmon_quant.sh:
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

