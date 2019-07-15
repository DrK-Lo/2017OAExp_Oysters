# Comparison between GTF files create from gff files directly from NCBI or those created by Kevin Johnson during FROGER workshop.

**Description**: Used `gffcompare NCBI_genome.GTF KM_genome.GTF` to compare the `.gtf` files that were created from the `.gff` file, both with the command `gffread`.

**Output**: Despite slight differences in the `.gff` files, the `.gtf` files are identical. This is due to the differences in the `.gff` files appead to be due to sequence types that are not retained in the `.gtf` files.  
