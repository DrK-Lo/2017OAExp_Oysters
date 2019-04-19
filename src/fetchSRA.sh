
#!/bin/bash

# Loop for downloading PE sequences from Accessions in my .txt file

for f in $F/*paired.text
do
    prefetch --option-file $f 
    while read -r LINE; do
        fastq-dump -O $F --split-files --readids $LINE
    done < $f
done
echo "STOP $date"