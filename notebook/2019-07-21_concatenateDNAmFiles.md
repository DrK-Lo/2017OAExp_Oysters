# Concatenate raw DNAm Reads

The reads from genewiz were separated across three lanes for each individual. These needed to be concantenated into a single file for each samples and read. In addition, these samples were with a general label, `{zr2576_1...zr2576_24}, rather than the sample labels specific to the 2017 exposure experiment (i.e. 17005), so they were renamed to the project specific labels using a [key](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/input_files/DNAm/DNA_seqKey_rv.csv) based on the sample sheet from the Roberts lab.

Input Files:
`zr2576_10_s1_R1.fastq.gz`
`zr2576_10_s2_R1.fastq.gz`
`zr2576_10_s3_R1.fastq.gz`

Example output Files:
`17005_DNAm_R1.fastq.gz`

Sample Function (concatenates single sample for read 1):
```
cat /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/${i}_*R1.fastq.gz > /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/${j}_DNAm_R1.fastq.gz
```

File Directory: `/shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp`

[Full Bash Script](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/src/cluster_scripts/DNAm/concate_rename.sh)
