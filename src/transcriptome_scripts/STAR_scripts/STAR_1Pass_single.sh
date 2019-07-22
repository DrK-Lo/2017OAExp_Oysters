#!/bin/bash

echo "Please put in raw file directory:"
read raw
echo "Sample ID (e.g. 17005):"
read sampleID
echo "Please put in name of new folder for output"
read output

base="/shared_lab/20180226_RNAseq_2017OAExp/RNA/"

echo "Outputs saving to : " $base$output

if [ -d "$base$output" ]; then
    echo "Directory Already Exists, please use another name"
else
    echo "Directory Created"
    mkdir "$base$output"
fi

echo "Processing the following sample: "
echo ls $base$raw/*$sampleID.R1.*

# This will loop through each sample in the raw folder directory
i=$base$raw/*$sampleID.R1.*
j=$base$raw/*$sampleID.R2.*

echo $i and $j
                      
                        /shared_lab/scripts/STAR --runThreadN 19 \
                        --genomeDir /shared_lab/20180226_RNAseq_2017OAExp/RNA/references/star_ref \
                        --outFilterMatchNminOverLread 0.17 --outFilterScoreMinOverLread 0.17 \
                        --readFilesIn $i $j \
                        --outSAMmapqUnique 40 \
                        --outSAMtype BAM Unsorted SortedByCoordinate \
                        --outFileNamePrefix $base$output/$sampleID_m2_ \
                        --readFilesCommand zcat
