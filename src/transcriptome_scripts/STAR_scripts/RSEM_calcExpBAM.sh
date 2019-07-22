#!/bin/bash

rsem-calculate-expression --alignments --paired-end --output-genome-bam -p 5 \
/shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_files/testBAM/m3/r2Aligned.toTranscriptome.out.bam \
/shared_lab/20180226_RNAseq_2017OAExp/RNA/references/RSEM_ref2/RSEM_ref2 \
/shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_files/testBAM/17005_

