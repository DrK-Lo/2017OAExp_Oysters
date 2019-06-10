# Revised STAR run folder `run20190610`

This is the latest STAR run that contains both gene counts outputs using `--quantiMode` in STAR and the modified parameters originally set by Brett.

**Input Info**

Raw Reads: `/shared*/2018*/RNA/rawfiles/*.R1.*`
* Example file: `P1_17162.R1.fq.gz`

Reference Genome Folder (created specifically for STAR): `/shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/CV_genome/star_wGTF`
* Note: this was created from the GCF.gtf file in the `CV_genome` folder and is a `.gtf` file created from `GCF_002022765.2_C_virginica-3.0_genomic.gff` using gffread.

Scripts: `/shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/scripts`

**Output Info**

Creates sequence junctions (sj), transcript read info (BAM files), and gene count data (readCount.out.tab) located in the `run20190610` folder.
* Note: `_m2_` indicates first pass outputs and `_m3_` for second pass (use second pass outputs for all downstream analysis.

### First Pass

Command line for performing first pass by running `bash_WholeSampleMapStar1Pass.sh`:
```
STAR_genomeMap_Whole_1Pass
```
* Note: this currently **isn't** a particularly flexible script and will need to be modified for each run. 
* This script should perform first pass on each sample within rawRead directory.

Bash Code for `bash_WholeSampleMapStar1Pass.sh`:
```
#!/bin/bash

for i in $( ls /shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/*.R1.* ); do
        for j in $( ls /shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/*.R2.* ); do
                file1=$( echo $i | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1)
                file2=$( echo $j | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1)

                if [ "$file1" == "$file2" ]
                then
                    	echo $i and $j
                        echo RNA"$file1"_m2
                        echo yes these files match
                        /shared_lab/scripts/STAR --runThreadN 19 \
                        --genomeDir /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/CV_genome/star_wGTF \
                        --outFilterMatchNminOverLread 0.17 --outFilterScoreMinOverLread 0.17 \
                        --readFilesIn $i $j --outSAMmapqUnique 40 --sjdbGTFtagExonParentGene gene_name \
                        --outSAMtype BAM Unsorted SortedByCoordinate \
                        --outFileNamePrefix /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/output/run20190610/"$file1"_m$
                        --readFilesCommand zcat
                fi
        done
done
```

### Second Pass

Command line for performing first pass by running `bash_WholeSampleMapStar2Pass.sh`:
```
STAR_genomeMap_Whole_2Pass
```
* Note: this currently **isn't** a particularly flexible script and will need to be modified for each run. 
* This script should perform first pass on each sample within rawRead directory.

Bash Script for :
```
#!/bin/bash

for i in $( ls /shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/*.R1.* ); do
        for j in $( ls /shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/*.R2.* ); do
                file1=$( echo $i | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1)
                file2=$( echo $j | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1)
                if [ "$file1" == "$file2" ]
                then
                    	m2_files=$( ls /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/output/run20190610/m2/*m2_SJ.out.tab)
                        echo $i and $j
                        echo RNA"$file1"_m3
                        echo yes these files match
                        echo $m2_files
                        /shared_lab/scripts/STAR --runThreadN 19 \
                        --genomeDir /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/CV_genome/star_wGTF \
                        --readFilesIn $i $j \
                        --outSAMmapqUnique 40 \
                        --outSAMtype BAM Unsorted SortedByCoordinate \
                        --sjdbGTFtagExonParentGene gene_name \
                        --quantMode GeneCounts --limitSjdbInsertNsj 1300000 \
                        --outFileNamePrefix /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/output/run20190610/"$file1"_m3_ \
                        --readFilesCommand zcat \
                        --sjdbFileChrStartEnd $m2_files
                fi
        done
done
```

