#!/bin/bash
  
PATH=$1
  
if [ ! -z "$PATH" ]
then
    wget --recursive -e robots=off --reject "index.html" --no-host-directories --no-parent https://ftp.ncbi.nih.gov/genomes/Crassostrea_virginica/ 
      
    $PATH
      
    mv $PATH/2017OAExp_Oysters/markdown_files/genomes/Crassostrea_virginica /home/downeyam/Github/2017OAExp_Oysters/Genome 
      
    cd $PATH/2017OAExp_Oysters/Genome/Assembled_chromosomes/seq
      
    cd $PATH/2017OAExp_Oysters/Genome/Assembled_chromosomes/seq
    zcat 6565_ref_C_virginica-3.0_chrMT.fa.gz 6565_ref_C_virginica-3.0_chr[0-9]*.fa.gz | sed 's/ref|//g' | sed 's/| Crass/  Crass/g' > /home/downeyam/Github/2017OAExp_Oysters/genome.fasta
      
    cd $PATH/2017OAExp_Oysters
    samtools faidx genome.fasta
    mawk -v OFS='\t' {'print $1,$2'} genome.fasta.fai > genome.file
      
    cd $PATH/2017OAExp_Oysters/Genome/GFF
    ln -s $PATH/2017OAExp_Oysters/Genome/Assembled_chromosomes/seq/genome.file .
else
    echo "ERROR:"
    echo "Add the path of the local directory following the file"
    exit 1
fi
