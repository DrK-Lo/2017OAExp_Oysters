#!/bin/bash

rsem-calculate-expression --star \
--paired-end \
--star-gzipped-read-file \
-p 8 \
/shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/P1_17005.R1.fq.gz \
/shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/P1_17005.R2.fq.gz \
/shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/CV_genome/RSEM_ref/RSEM_ref \
/shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/output/run20190610/RSEM

