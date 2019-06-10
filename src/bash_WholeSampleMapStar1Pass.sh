#!/bin/bash
for i in $( ls /shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/*.R1.* ); do
        for j in $( ls /shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/*.R2.* ); do
                file1=$( echo $i | cut -d'_' -f 1)
                file2=$( echo $j | cut -d'_' -f 1)
                if [ "$file1" == "$file2" ]
                then
                        echo $i and $j
                        echo RNA"$file1"_m2
                        echo yes these files match
                        /shared_lab/scripts/STAR --runThreadN 32 \
                        --genomeDir /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/CV_genome/star \
                        --readFilesIn $i $j \
                        --quantMode TranscriptomeSAM GeneCounts \
                        --outSAMtype BAM Unsorted SortedByCoordinate \
                        --outFileNamePrefix /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/output/rerun/RNA"$file1"_m2 \
                        --readFilesCommand zcat
                fi
        done
done