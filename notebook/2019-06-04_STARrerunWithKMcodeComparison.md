# Star alignment comparison (first pass) between original alignment and using code from Kevin Johnson.

**Overview**: I reran the STAR alignment based on a newly asembled reference genome from code borrowed from Kevin, then compared those alignments to the mapping that was originally perfomred by Brett.

### Creating new genome file

A new STAR genome file was created utilizing the latest `.fna` fasta format and the `.gff` file from the NCBI website to create the reference with code from Kevin Johnson.

Sample of original reference code for generating genome files: 

```
/home/kjohnson/STAR/bin/Linux_x86_64/STAR --runThreadN 40 \
--runMode genomeGenerate \
--genomeDir /work/kjohnson/CV_genome/star \
--genomeFastaFiles /work/kjohnson/CV_genome/GCF_002022765.2_C_virginica-3.0_genomic.fna \
--sjdbGT /work/kjohnson/CV_genome/GCF_002022765.2_C_virginica-3.0_genomic.gff 
```

**NCBI FILES**: 
Fasta file used: GCF_002022765.2_C_virginica-3.0_genomic.fna
Gff file used: GCF_002022765.2_C_virginica-3.0_genomic.gff

**Notes**:
Everything else was default. Also as a warning getting this to run was a little challenging, and it kept throwing a `do not recognize the directory` message that was due to an additional space in the script after `--genomeFastaFiles`, be sure to remove any additional spaces before running.

### First Pass Mapping 

With the genome files made I reran the first pass of the STAR mapper on sample **17005**.

Sample of original code:

```
/home/kjohnson/STAR/bin/Linux_x86_64/STAR --genomeDir /work/kjohnson/CV_genome/star \
--readFilesCommand zcat \
--outSAMmapqUnique 60 \
--readFilesIn /work/kjohnson/CVirginica_transcriptome/Trim/5_NB_30_F_paired.fq.gz /work/kjohnson/CVirginica_transcriptome/Trim/5_NB_30_R_paired.fq.gz \
--runThreadN 16 \
--outFileNamePrefix “/work/kjohnson/CVirginica_transcriptome/Star/5_NB_30”
```

**Read files used in actual script** (found in the `RNA/rawfiles` folder):
* `/P1_17005.R1.fq.gz`  
* `/P1_17005.R1.fq.gz`  

## File Location.
All files are stored on `comp5` in the `/shared*/2018*/RNA/STAR_KM*/` directory, including the specific bash scripts I used to create the genome file and perform the mapping.

* `/output` : contians the alignment outputs
* `/scripts` : contians the bash scripts used based on Kevin Johnsons code.
* `/CV_genome`: contians `.fna`, `.gff`, and the assemblied STAR genome folders

## Comparison of STAR rerun with Kevin script verse orignal

![]()

Figure Description: Original run log for samples **17005** after the 1st pass in on the left, and the 1st pass using the new scripts is on the right.

**Thoughts**: Looks like we still had better alignment during the original run based on most log metrics. Importantly this include a higher % mapping for unique reads. This might indicate that the original run is still our best data, but ill be confirming this by comparing the actual outputs to see if these are identifying the same reads for each transcript feature.

**Additional Notes**: With this test I was also able to confirm that the `P1_17***.R1.fq.gz` files are indeed the post dDocent and trimmomatic files that Brett used for his STAR mapping run based on the mathcing number of input reads in both of our STAR runs.
