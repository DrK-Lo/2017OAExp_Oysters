#!/bin/bash

salmon quant -i oyster_index_from_genome \
	-l A \
	-1 /shared_lab/20180226_RNAseq_2017OAExp/raw/17005_R1_001.fastq.gz \
	-2 /shared_lab/20180226_RNAseq_2017OAExp/raw/17005_R2_001.fastq.gz \
	-p 10 \
	--validateMappings \
	--rangeFactorizationBins 4 \
	--numBootstraps 10 \
	--gcBias \
	--seqBias \
	-o test_out/17005_fromGenome_bias2
