#!/bin/bash

#echo "Please put in raw file directory:"
#read raw
#echo "Sample ID (e.g. 17005):"
#read sampleID
#echo "Please put in name of new folder for output"
#read output

base="/shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles"

#echo "Outputs saving to : " $base$output

#if [ -d "$base$output" ]; then
#    echo "Directory Already Exists, please use another name"
#else
#    echo "Directory Created"
#    mkdir "$base$output"
#fi

#echo "Processing the following sample: "
#echo ls $base$raw/*$sampleID.R1.*

# This will loop through each sample in the raw folder directory
i=$base/*17005.R1.*
j=$base/*17005.R2.*

echo $i and $j
m2_files=$(ls /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_files/testBAM/m2/SJ.out.tab)
                      
                        /shared_lab/scripts/STAR --runThreadN 19 \
                        --genomeDir /shared_lab/20180226_RNAseq_2017OAExp/RNA/references/star_ref \
                        --readFilesIn $i $j \
                        --outSAMmapqUnique 40 \
                        --quantMode TranscriptomeSAM GeneCounts --limitSjdbInsertNsj 1500000 \
                        --outSAMtype BAM Unsorted SortedByCoordinate \
                        --outFileNamePrefix /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_files/testBAM/m3/r2 \
                        --readFilesCommand zcat

