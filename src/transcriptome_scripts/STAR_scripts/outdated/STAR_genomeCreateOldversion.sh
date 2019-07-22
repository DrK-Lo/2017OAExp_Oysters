!#/bin/bash

STAR --runThreadN 32 \
--runMode genomeGenerate \
--genomeDir /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/CV_genome/star_wGTF \
--genomeFastaFiles /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/CV_genome/GCF_002022765.2_C_virginica-3.0_genomic.fna \
--sjdbGTFfile /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/CV_genome/GCF.gtf
