#!/bin/bash

files = ls /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp

echo Files being trimmed: $files

trim_galore \
--paired \
--clip_r1 10 \
--clip_r2 10 \
--three_prime_clip_R1 10 \
--three_prime_clip_R2 10 \
--output_dir /shared_lab/20180226_RNAseq_2017OAExp/DNAm/20190719_fastqc_trim_10bp_Cvirginica_MBD \
--fastqc_args "--outdir /shared_lab/20180226_RNAseq_2017OAExp/DNAm/20190719_trim_galore_files --threads 18" \
$files \
2> /shared_lab/20180226_RNAseq_2017OAExp/DNAm/20190719_fastqc_trim_10bp_Cvirginica_MBD/stderr.log
