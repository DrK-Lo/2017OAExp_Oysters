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
