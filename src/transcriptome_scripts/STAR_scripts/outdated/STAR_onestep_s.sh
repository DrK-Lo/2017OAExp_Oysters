!#/bin/bash

STAR --genomeDir /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/CV_genome/star_wGTF \
--readFilesCommand zcat --outSAMmapqUnique 60 \
--readFilesIn /shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/P1_17005.R1.fq.gz /shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/P1_17005.R2.fq.gz \
--runThreadN 15 --twopassMode Basic \
--sjdbGTFtagExonParentGene gene_name \
--quantMode GeneCounts \
--outFileNamePrefix /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/output/17005_GeneCounts/17005_
