#!/bin/bash

salmon quant -i oyster_index_from_genome \
	-l A \
	-1 /shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/P1_17005.F.fq.gz \
	-2 /shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/P1_17005.R.fq.gz \
	-p 10 \
	--validateMappings \
	--rangeFactorizationBins 4 \
	--numBootstraps 10 \
	-o test_out/17005_fromRawFiles_F2_fromGenomeIndex
