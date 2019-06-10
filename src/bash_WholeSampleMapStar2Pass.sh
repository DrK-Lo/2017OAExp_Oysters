#!/bin/bash  

for i in $( ls /shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/*.R1.* ); do
        for j in $( ls /shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/*.R2.* ); do
                file1=$( echo $i | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1)
                file2=$( echo $j | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1)
                if [ "$file1" == "$file2" ]
                then
                        m2_files=$( ls /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/output/rerun3/*m2_SJ.out.tab)
                        echo $i and $j
                        echo RNA"$file1"_m3
                        echo yes these files match
                        echo $m2_files
                        /shared_lab/scripts/STAR --runThreadN 32 \
                        --genomeDir ./genome/ \
                        --readFilesIn $i $j \
                        --outSAMmapqUnique 40 \
                        --outSAMtype BAM Unsorted SortedByCoordinate \
                        --sjdbGTFtagExonParentGene gene_name \
                        --quantMode GeneCounts \
                        --outFileNamePrefix /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_KM_version/output/rerun3/"$file1"_m3_ \
                        --readFilesCommand zcat \
                        --sjdbFileChrStartEnd "$m2_files"
                fi
        done
done