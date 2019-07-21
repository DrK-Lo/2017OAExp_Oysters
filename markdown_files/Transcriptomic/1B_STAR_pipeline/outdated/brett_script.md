Cvirginica\_palialpH\_RNAseq\_bioinformatics
================
Brett Ford
5/22/2018

## From Raw Sequence Data to Count Matrix

This document details the computational steps needed to go from raw
RNAseq sequences from the Eastern oyster (*Crassostrea virginica*) to a
count matrix, with the ultimate goal of identifying differential gene
expression between treatment groups.

The general outline consists of assembling the genome, trimming and
removing adapters from RNAseq reads, mapping reads to the assembled
genome, quality filtering, and then creating a count matrix from high
quality mapped reads.

The following code assumes that you have created a project directory
with the subdirectory ./RNA (with raw RNA sequences)

### Assemble Genome

Used the following code to download the genome to your project directory
(replace last folder extension name project
folder)

``` bash
wget --recursive -e robots=off --reject "index.html" --no-host-directories --no-parent
https://ftp.ncbi.nih.gov/genomes/Crassostrea_virginica/ -P 
/shared_lab/20180226_RNAseq_2017OAExp
```

This will likely create a directory called
`/shared_lab/20180226_RNAseq_2017OAExp/genomes` with the subdirectory
`/shared_lab/20180226_RNAseq_2017OAExp/genomes/Crassostrea_virginica`.
Move the Crassostrea\_viriginica directory to the project directory and
rename to Genome to follow Jon???s code

``` bash
mv /shared_lab/20180226_RNAseq_2017OAExp/genomes/Crassostrea_virginica
/shared_lab/20180226_RNAseq_2017OAExp
cd /shared_lab/20180226_RNAseq_2017OAExp
mv Crassostrea_virginica Genome
rm -rf genomes
```

Change to the directory containing sequencing information for all
chromosomes

``` bash
cd /shared_lab/20180226_RNAseq_2017OAExp/Genome/Assembled_chromosomes/seq
```

Assemble full reference genome from all
chromosomes

``` bash
zcat 6565_ref_C_virginica-3.0_chrMT.fa.gz 6565_ref_C_virginica-3.0_chr[0-9]*.fa.gz | 
sed 's/ref|//g' | sed 's/| Crass/  Crass/g' > genome.fasta
```

Create a genome file to use with bedtools sorting

``` bash
samtools faidx genome.fasta
mawk -v OFS='\t' {'print $1,$2'} genome.fasta.fai > genome.file
```

Change to the genome features file (GFF) directory and set up link to
genome file

``` bash
cd /shared_lab/20180226_RNAseq_2017OAExp/Genome/GFF
ln -s /shared_lab/20180226_RNAseq_2017OAExp/Genome/Assembled_chromosomes/seq/genome.file .
```

### Trim Reads and Remove Adapters

Change to subdirectory with the raw sequence information and download
the adapter file corresponding to your library prep kit.

``` bash
cd /shared_lab/20180226_RNAseq_2017OAExp/RNA
wget https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE-2.fa
```

Download custom dDocent v2.2.ed20 that looks for adapters in working
directory

``` bash
wget https://raw.githubusercontent.com/jpuritz/EecSeq/master/Bioinformatics/dDocent
```

Download the configuration
file

``` bash
wget https://raw.githubusercontent.com/jpuritz/EecSeq/master/Bioinformatics/RNA.config
```

Use nano to edit config file to appropriate settings (Number of
Processors=45, Maximum Memory=100G, Trimming=yes, Assembly?=no,
Type\_of\_Assembly= , Clustering\_Similarity%= , Mapping Reads?=no,
Mapping\_Match\_Value= , Mapping\_Mismatch\_Value= ,
Mapping\_GapOpen\_Penalty= , Calling\_SNPs?=no,
<Email=br.ford@northeastern.edu>)

``` bash
nano RNA.config
```

Use nano to edit dDocent to trim accoding to custom adapter file (use
ctrl+w to find and replace TruSeq2-PE.fa with TruSeq3-PE-2.fa; 2
instances)

``` bash
nano dDocent
```

Make sure that the input files are in the format specified by dDocent.
That is, filenames should have the extensions \*.R1.fq.gz for forward
reads and \*.R2.fq.fz for reverse reads (see <http://ddocent.com/quick/>
for more details).

Run dDocent. You may have to change permissions to run it.

``` bash
./dDocent RNA.config
```

#### What’s going on behind the code?

The dDocent pipeline uses Trimmomatic trimming tool to remove adapter
and low quality sequences from the ends of reads. Within the dDocent
code, it is specified to be paired-end (which is automatically
recognized based on our file naming scheme), removes adapters based on
thresholds for how well the adapter sequences align to reads (2:30:10;
see Trimmomatic manual for more details), removes leading bases with
phred quality score less than 20, removes trailing bases with phred
quality score less than 20, scans the reads at a 5bp window and cuts
when the average quality of the five bases is less than 10, and makes
sure all reads are a minimum length after this cutting (greater than the
shortest read/2).

### Map Reads to Genome

Make sure you are still in the RNA subdirectory. Create a link to the
reference genome (the code below is different from original code from
EecSeq bioinformatics pipeline because of naming above) and a link to
the gff file. Also make sure to create a `./genome` directory. Make sure
STAR is installed on your system. You may need to specify path to STAR.
Make sure you get ouput when calling STAR, otherwise you need to install
it on your system (see here: <https://github.com/alexdobin/STAR> for
installation instructions)

``` bash
cd /shared_lab/20180226_RNAseq_2017OAExp/RNA
ln -s /shared_lab/20180226_RNAseq_2017OAExp/Genome/Assembled_chromosomes/seq/genome.fasta
reference.fasta
ln -s /shared_lab/20180226_RNAseq_2017OAExp/Genome/GFF/ref_C_virginica-3.0_top_level.gff3 
top.gff
mkdir genome
/shared_lab/scripts/STAR
```

Create genome files (genome index) with STAR
aligner

``` bash
STAR --runMode genomeGenerate --runThreadN 64 --genomeDir ./genome --genomeFastaFiles
reference.fasta --sjdbGTFfile top.gff --sjdbGTFtagExonParentTranscript Parent 
--sjdbOverhang 149
```

The code above is fairly self explanatory. We are creating a genome
index from the reference genome in our genome directory, specifying our
annotations, which are in gff format. The
`--sjdbGTFtagExonParentTranscript Parent` argument is required for gff3
annotation formats (the format of our top.gff file). The
`--sjdbOverhang 149` specifies that the genomic sequence must be 149 bp
in length around the annotated junction to be used in constructing the
splice junctions database. Ideally, the length is readlength-1, hence
150-1=149.

To map our reads, we will use the STAR aligner. STAR uses a 2-pass
mapping system where we first map each individual sample to the genome,
and then we conduct a second mapping pass for all samples (see
'Multi-sample 2-pass mapping' in the STAR manual)

Run STAR in dual-pass mode. You can either align reads from each
individual to the genome
using:

``` bash
STAR --runThreadN 32 --genomeDir ./genome/ --readFilesIn  17005_R1_001.fastq.gz
17005_R2_001.fastq.gz  --outFilterMatchNminOverLread 0.17 --outFilterScoreMinOverLread 
0.17 --outSAMmapqUnique 40 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix 
RNA17005_m2  --readFilesCommand zcat
```

…and then repeat for each individual, OR use for loop to loop through
all fastq files in directory using:

``` bash
#!/bin/bash
for i in $( ls *_R1_* ); do
        for j in $( ls *_R2_* ); do
                file1=$( echo $i | cut -d'_' -f 1)
                file2=$( echo $j | cut -d'_' -f 1)
                if [ "$file1" == "$file2" ]
                then
                        echo $i and $j
                        echo RNA"$file1"_m2
                        echo yes these files match
                        /shared_lab/scripts/STAR --runThreadN 32 \
                        --genomeDir ./genome/ \
                        --readFilesIn $i $j --limitSjdbInsertNsj 1620452 \
                        --outFilterMatchNminOverLread 0.17 \
                        --outFilterScoreMinOverLread 0.17 \
                        --outSAMmapqUnique 40 --outSAMtype BAM SortedByCoordinate \
                        --outFileNamePrefix RNA"$file1"_m2  --readFilesCommand zcat
                fi
        done
done
```

The `--limitSjdbInsertNsj` argument specifies the maximum number of
junction to be inserted to the genome on the fly at the mapping stage. I
had to increase this from the default of 1000000 for the aligner to
properly work. The `--outFilterMatchNminOverLread` and
`--outFilterScoreMinOverLread` specifies that the alignment will be
output if the number of matched bases is higher than or equal to this
value (normalized to the read length, which in this case is sum of
mates’ lengths for paired-end reads; 0.17 allows reads that align at
least 51 (0.17 is 51% of 300) of their bases to the 300 bp of mates’
lengths) and whether the score is higher than or equal to this value,
respectively. These values should be the same. The
`--outSAMmapqUnique 40` argument specifies that the MAPQ value for the
sam output should have a max of 40 (there are different threshold
depending on the aligners).

Now run the second pass, running the following
code:

``` bash
/shared_lab/scripts/STAR --runThreadN 32 --genomeDir ./genome/ --readFilesIn  17005_R1_001.fastq.gz 17005_R2_001.fastq.gz  --limitSjdbInsertNsj 1620452 --outFilterMatchNminOverLread 0.17 --outFilterScoreMinOverLread 0.17 --outSAMmapqUnique 40 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix RNA17005_m3  --readFilesCommand zcat --sjdbFileChrStartEnd ./RNA17005_m2SJ.out.tab ./RNA17007_m2SJ.out.tab ./RNA17013_m2SJ.out.tab ./RNA17019_m2SJ.out.tab ./RNA17069_m2SJ.out.tab ./RNA17070_m2SJ.out.tab ./RNA17072_m2SJ.out.tab ./RNA17079_m2SJ.out.tab ./RNA17090_m2SJ.out.tab ./RNA17094_m2SJ.out.tab ./RNA17099_m2SJ.out.tab ./RNA17108_m2SJ.out.tab ./RNA17122_m2SJ.out.tab ./RNA17130_m2SJ.out.tab ./RNA17142_m2SJ.out.tab ./RNA17145_m2SJ.out.tab ./RNA17162_m2SJ.out.tab ./RNA17174_m2SJ.out.tab ./RNA17176_m2SJ.out.tab ./RNA17178_m2SJ.out.tab ./RNA17181_m2SJ.out.tab ./RNA17203_m2SJ.out.tab ./RNA17211_m2SJ.out.tab ./RNA17213_m2SJ.out.tab
```

… changing the `--readFilesIn` argument to each individual, OR by using
a for loop once again:

``` bash
#!/bin/bash
for i in $( ls *_R1_* ); do
        for j in $( ls *_R2_* ); do
                file1=$( echo $i | cut -d'_' -f 1)
                file2=$( echo $j | cut -d'_' -f 1)
                if [ "$file1" == "$file2" ]
                then
                        m2_files=$( ls *m2SJ.out.tab)
                        echo $i and $j
                        echo RNA"$file1"_m3
                        echo yes these files match
                        echo $m2_files
                        /shared_lab/scripts/STAR --runThreadN 32 \
                        --genomeDir ./genome/ \
                        --readFilesIn $i $j --limitSjdbInsertNsj 1620452 \
                        --outFilterMatchNminOverLread 0.17 \
                        --outFilterScoreMinOverLread 0.17 \
                        --outSAMmapqUnique 40 --outSAMtype BAM SortedByCoordinate \
                        --outFileNamePrefix RNA"$file1"_m3  --readFilesCommand zcat \
                        --sjdbFileChrStartEnd "$m2_files"
                fi
        done
done
```

The first pass of STAR finds splice junction loci and the second pass
aligns spliced reads with short overhangs across the junctions that were
identified in the first pass.

### Filter Reads that Poorly Mapped to Genome

To filter reads that poorly mapped, you can either run the following
code each
individual:

``` bash
samtools view -@64 -q4 -h -F 0x100 -F 0x400 RNA17005_m3Aligned.sortedByCoord.out.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/' | samtools view -b > RNA17005.q4.bam 2> RNA17005.q4_error.log
```

… replacing the individual name, OR use the following bash script:

``` bash
#!/bin/bash
for i in $( ls RNA*m3Aligned.sortedByCoord.out.bam ); do
        ind=$( echo $i | cut -d'_' -f 1)
        samtools view -@64 -q4 -h -F 0x100 -F 0x400 "$i"| \
        mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/' | \
        samtools view -b > "$ind".q4.bam 2> "$ind".q4_error.log
done
```

The argument -q4 skips alignments with MAPQ smaller than 4, the -h
argument specifies to include header, -F specifies to not output
alignments with with any flag between 100 and 400, and the mawk line
filters out reads that had significant hard or soft clipping \>79 bp.

\#\#\#Post Processing- Creating Count Matrix First, sort all bam files
by name to use in htseq-count. Sorting by name seems to work the best
when working with bam files output from STAR alignments. You’ll need to
create a tmp file to hold temporary files required for samtools

``` bash
mkdir tmp
```

Then sort files:

``` bash
#!/bin/bash
for i in $( ls RNA*.q4.bam ); do
        ind=$( echo $i | cut -d'.' -f 1)
        samtools sort -n -@32 -T /tmp/"$i".q4.namesorted \
        -o "$ind".q4.namesorted.bam $i 2> "$ind"_bam_namesorted_error.log
done
```

I then converted the gff file to a gtf file using gffread (see here:
<https://github.com/gpertea/gffread> for installation instructions),
because it was easier to use gtf file for htseq, the program that gives
a count matrix. Make sure that gffread executable is in the current
directory (`./RNA` subdirectory).

``` bash
./gffread ref_C_virginica-3.0_top_level.gff3 -T -o top.gtf
```

At this point there may be missing gene\_name data in the `top.gtf`
file. Because the gene\_name is the id attribute used to count
transcripts in the code below, it is essential that each element has a
gene\_name. Therefore, I manually added a gene\_name for elements
without one starting “unknown\_gene1” all the way up to
“unknown\_gene26”

``` bash
nano top.gtf
```

Now use the name sorted bam files and the gtf file in htseq-count to
count the average number of transcripts (exons specifically) per gene.

``` bash
#!/bin/bash
for i in $( ls RNA*.q4.namesorted.bam ); do
        ind=$( echo $i | cut -d'.' -f 1| cut -c 4-)
        htseq-count -f bam -r name -s yes -i gene_name \
        --additional-attr gene_id -t exon -m union \
        $i top.gtf > \
        "$ind"_exon_count.txt 2> "$ind"_exon_count_error.log
done
```

This creates files for each individual. I then created an Rscript to
combine all of these files into a single count matrix with each
individual as a separate column and each gene as a separate row.

``` r
base_path <- ("/shared_lab/20180226_RNAseq_2017OAExp/RNA")
exon_count_filenames <- list.files(path = base_path, pattern = ".+_exon_count.txt")
k <- 1
column_names <- c("gene_name", "gene_id")
for (i in exon_count_filenames) {
    file_contents <- unlist(strsplit(i, split = "[_]"))
    ind <- file_contents[1]
    table_name <- paste0(ind, "_exon_count")
    assign(table_name, read.table(i, header = F, sep = "\t"))
    if (k == 1) {
        C_virginica_gene_count <- data.frame(get(table_name))
    } else {
        C_virginica_gene_count <- cbind(C_virginica_gene_count, get(table_name)$V3)
    }
    k = k + 1
    column_names <- c(column_names, paste0("RNA", ind))
}
names(C_virginica_gene_count) <- column_names
write.table(C_virginica_gene_count, file = "C_virginica_gene_count.txt", col.names = T)
```

\#\#Additional code to visualize overall mapping efficiency

#### Create bed files to compare coverage after mapping reads

Make sure BEDOPS is installed on your system (see here:
<https://github.com/bedops> for installation instructions). Convert gff
file to bed file and then reduce to a sorted bed file with only
exons.

``` bash
gff2bed < ref_C_virginica-3.0_top_level.gff3 > ref_C_virginica-3.0_top.bed
mawk '$8 ~ /exon/' ref_C_virginica-3.0_top.bed > ref3.0.exon.bed
bedtools sort -i ref3.0.exon.bed -faidx <(cut -f1 genome.file) > sorted.ref3.0.exon.bed
```

Remove duplicates and concatenate overlapping intervals from transcript
variants

``` bash
bedtools merge -i sorted.ref3.0.exon.bed > sorted.ref3.0.exon.sc.bed
```

Create bed file for all gene regions

``` bash
mawk '$8 ~ /gene/' ref_C_virginica-3.0_top.bed > ref3.0.gene.bed
bedtools sort -i ref3.0.gene.bed -faidx <(cut -f1 genome.file) > sorted.ref3.0.gene.bed
bedtools merge -i sorted.ref3.0.gene.bed > sorted.ref3.0.gene.sc.bed
```

Create bed file for intergenic, intron, non-coding, and CDS
regions

``` bash
bedtools complement -i sorted.ref3.0.gene.bed -g genome.file -sorted | bedtools subtract -b sorted.ref3.0.exon.sc.bed -a - > cv.ref3.intergenic.bed
bedtools complement -i sorted.ref3.0.exon.sc.bed -g genome.file -sorted > cv.ref3.noncoding.bed
bedtools intersect -a cv.ref3.noncoding.bed -b sorted.ref3.0.gene.sc.bed -sorted > cv.ref3.intron.bed
mawk '$8 ~ /CDS/' ref_C_virginica-3.0_top.bed > ref3.0.CDS.bed
bedtools sort -i ref3.0.CDS.bed -faidx <(cut -f1 genome.file) > sorted.ref3.0.CDS.bed
bedtools merge -i sorted.ref3.0.CDS.bed > sorted.ref3.0.CDS.sc.b
bedtools sort -i sorted.ref3.0.CDS.sc.b -faidx <(cut -f1 genome.file) > sorted.ref3.0.CDS.sc.bed
```

Create bed for untranslated regions (UTRs) of
exons

``` bash
bedtools subtract -a sorted.ref3.0.exon.bed -b sorted.ref3.0.CDS.bed -g genome.file -sorted > sorted.ref3.0.UTR.bed
bedtools merge -i <(bedtools sort -i sorted.ref3.0.UTR.bed -faidx <(cut -f1 genome.file))> sorted.ref3.0.UTR.sc.b
bedtools sort -i sorted.ref3.0.UTR.sc.b -faidx <(cut -f1 genome.file) > sorted.ref3.0.UTR.sc.bed
```

Create bed file for mtDNA

``` bash
mawk '$1 ~ /NC_007175.2/' ref_C_virginica-3.0_top.bed > mtDNA.bed
```

Merge all alignment files into a single merged
file

``` bash
samtools merge -@64 m4.merged.bam RNA17005_m3Aligned.sortedByCoord.out.bam RNA17007_m3Aligned.sortedByCoord.out.bam RNA17013_m3Aligned.sortedByCoord.out.bam RNA17019_m3Aligned.sortedByCoord.out.bam RNA17069_m3Aligned.sortedByCoord.out.bam RNA17070_m3Aligned.sortedByCoord.out.bam RNA17072_m3Aligned.sortedByCoord.out.bam RNA17079_m3Aligned.sortedByCoord.out.bam RNA17090_m3Aligned.sortedByCoord.out.bam RNA17094_m3Aligned.sortedByCoord.out.bam RNA17099_m3Aligned.sortedByCoord.out.bam RNA17108_m3Aligned.sortedByCoord.out.bam RNA17122_m3Aligned.sortedByCoord.out.bam RNA17130_m3Aligned.sortedByCoord.out.bam RNA17142_m3Aligned.sortedByCoord.out.bam RNA17145_m3Aligned.sortedByCoord.out.bam RNA17162_m3Aligned.sortedByCoord.out.bam RNA17174_m3Aligned.sortedByCoord.out.bam RNA17176_m3Aligned.sortedByCoord.out.bam RNA17178_m3Aligned.sortedByCoord.out.bam RNA17181_m3Aligned.sortedByCoord.out.bam RNA17203_m3Aligned.sortedByCoord.out.bam RNA17211_m3Aligned.sortedByCoord.out.bam RNA17213_m3Aligned.sortedByCoord.out.bam
```

Again, filter reads that did not uniquely map to the genome from this
merged
file

``` bash
samtools view -@64 -q4 -h -F 0x100 -F 0x400 m4.merged.bam| mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/' | samtools view -b > m4.q4.merged.bam
```

Now we want to look at coverage of our reads across the genome. To do
so, first create links to the bed files in the genome features file
(GFF) we downloaded earlier to see can look at coverage for various
genome regions.

``` bash
ln -s /shared_lab/20180226_RNAseq_2017OAExp/Genome/GFF/*.bed .
ln -s /shared_lab/20180226_RNAseq_2017OAExp/Genome/GFF/genome.file .
```

Then, let’s get some basic numbers on the number of reads mapping, the
number mapping to genes, and the number mapping to
exons:

``` bash
paste <(samtools view -@32 -c m4.q4.merged.bam) <(samtools view -@32 -c -L sorted.ref3.0.gene.bed m4.q4.merged.bam) <(samtools view -@32 -c -L sorted.ref3.0.exon.sc.bed m4.q4.merged.bam)
```

Next we’ll calculate per base pair coverage across genomic regions and
output the data in txt format:

``` bash
#introns
bedtools coverage -hist -b m4.q4.merged.bam -a cv.ref3.intron.bed -g genome.file -sorted  -split | grep ^all > AllRNAm4q4.hist.AllIntron.all.split.txt
#intergenic
bedtools coverage -hist -b m4.q4.merged.bam -a cv.ref3.intergenic.bed -g genome.file -sorted -split | grep ^all > AllRNAm4q4.hist.AllIntergenic.all.split.txt
#exons (genes-introns)
bedtools coverage -hist -b m4.q4.merged.bam -a sorted.ref3.0.exon.sc.bed -g genome.file -sorted  -split | grep ^all > AllRNAm4q4.hist.AllExon.all.split.txt
#untraslated regions
bedtools coverage -hist -b m4.q4.merged.bam -a sorted.ref3.0.UTR.sc.bed -g genome.file -sorted  -split | grep ^all > AllRNAm4q4.hist.allUTR.all.split.txt
#coding sequences (genes-introns-UTR)
bedtools coverage -hist -b m4.q4.merged.bam -a sorted.ref3.0.CDS.sc.bed -g genome.file -sorted  -split | grep ^all > AllRNAm4q4.hist.AllCDS.all.split.txt
```

Next we’ll use an Rscript to plot the coverage across genome regions.
This was modified from the original code on the EecSeq GitHub page:

``` r
setwd("/shared_lab/20180226_RNAseq_2017OAExp/RNA")
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)
print(files <- list.files(pattern = "AllRNAm4q4.hist."))
labs <- c("CDS", "Exon", "Intergenic", "Intron", "UTR")
cbPalette <- c("#D55E00", "#009E73", "#56B4E9", "#0072B2", "#E69F00", "#F0E442", 
    "#999999", "#CC79A7", "#7570B3")
cov <- list()
for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])[, c(2, 5)]
    cov_cumul = 1 - cumsum(cov[[i]][, 2])
    cov[[i]]$cov_cumul <- c(1, cov_cumul[-length(cov_cumul)])
    cov[[i]]$sample = labs[i]
}
cov_df = do.call("rbind", cov)
names(cov_df)[1:2] = c("depth", "fraction")
pcbPalette <- c("#009E73", "#D55E00", "#CC79A7", "#0072B2", "#56B4E9")
p <- ggplot(cov_df, aes(x = depth, y = cov_cumul, color = sample)) + xlim(0, 
    100) + scale_alpha(guide = "none") + geom_line(size = 1.5) + # geom_segment(aes(x=20, y=0, xend=20, yend=1, color='red'))+
scale_color_manual(values = pcbPalette) + scale_fill_manual(values = pcbPalette) + 
    ylab("% of Bases > Depth") + xlab("Depth") + theme_bw() + theme(plot.title = element_text(size = 12, 
    face = "bold", hjust = 0.5)) + theme(legend.title = element_blank()) + theme(legend.position = c(0.75, 
    0.75))
png(filename = "Figure2_BFedited.png", type = "cairo", units = "px", width = 5600, 
    height = 3000, res = 600, bg = "transparent")
p
dev.off()
```

![Figure1](/Users/brettford/Desktop/Northeastern/lab_work/20180119_RNAseq_proj/Figure2_BFedited.png)

Displays the percent of bases covered for each region type by your reads
We see that coding sequences, exons, and UTRs have the greatest
percentage of bases that are covered by reads, and that as depth
increases (depth here meaning the number of positions you are looking
at) increases, there are a lower percentage of bases covering them.
