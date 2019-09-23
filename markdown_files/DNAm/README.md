# DNA Methylation Analysis


## Creating cytosine count matrices from raw reads

### Step One - Downloading the Data

Data was sent to the Roberts lab after sequencing and stored on their server [(Link for details)](https://robertslab.github.io/sams-notebook/2019/06/26/Data-Received-C.virginica-Mantle-MBD-BSseq-from-ZymoResearch.html). Below is the code used to download the data from their server:
```
wget -r https://owl.fish.washington.edu/nightingales/C_virginica/
```
Then moved to the appropriate working directory: `/shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp`

### Step Two - Concatenating the Data
Paired end sequencing for each sample was run across three lanes. These needed to be first concatenated into a single file and renamed so that the file ID matched the sample ID used elsewhere. This was done in a single step, using the `cat` function to combine files and storing the output in a new file that was renamed using the sampleID.

Directory: `/shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp`  

Raw file example: `zr2576_10_s1_R1.fastq.gz`
  * `zr2576_10` : Sample ID from the sequencer
  * `s1` : sequence lane (i.e. `s1`-`s3`)
  * `R1` : read 1 (paired end reads, `R1` or `R2`)
  
You will also need to ensure the `DNA_seqKey._rv.csv` is up to date and links all raw file labels with the experimental sample labels.
Example: `zr2576_10      17094`
  
Bash script used to concatenate files:
```
#!/bin/bash

filename='/shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/DNA_seqKey_rv.csv'

while read line; do
        # reading each line
        echo "Matching DNAm ID with sample : $line"

        IFS='    ' # four spaces set as delimiter
        read -ra ADDR <<< "$line" # str is read into an array as tokens separated by IFS
        i="${ADDR[0]}" # Save first column
        j="${ADDR[1]}" # Save second column

        echo "Concatenating sample lists:"
        echo ls /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/${i}_*R1.fastq*
        echo "Output file:"
        echo "/shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/${j}_DNAm_R1.fastq.gz"

        cat /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/${i}_*R1.fastq.gz > /shared_lab/20180226_RNAseq_2017OAExp/DNAm/raw$
        cat /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/${i}_*R2.fastq.gz > /shared_lab/20180226_RNAseq_2017OAExp/DNAm/raw$
        echo ""
        echo ""
done < $filename
```
Script [Link]()

Outputs were saved in the same directory. Example of file output: `17203_DNAm_R1.fastq.gz`

### Step Three - Perform quality and adapter trim
Remove reads with poor mapping quality and also cut off 10bp from both the 5' and 3' regions of either strand. This will remove any adapter sequence to improve downstream mapping. This trimming follows the recommendations of bismark when the library prep was down with the pico zymo kit. 

Starts with the rawfiles created from the previous step.
Example: `/shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/17094_DNAm_R2.fastq.gz`

Bash script for quality trimming:
```
#!/bin/bash

files= ls /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/17*

echo "Files being trimmed: $files"


for i in $( ls /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/17*R1*); do
        for j in $( ls /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/17*R2*); do
                R1=$( echo $i | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 1)
                R2=$( echo $j | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 1)

                if [ "$R1" == "$R2" ]
                then
                    	echo "Starting with sample $R1"

                        trim_galore \
                        --paired \
                        --clip_r1 10 \
                        --clip_r2 10 \
                        --three_prime_clip_R1 10 \
                        --three_prime_clip_R2 10 \
                        --output_dir /shared_lab/20180226_RNAseq_2017OAExp/DNAm/20190719_fastqc_trim_10bp_Cvirginica_MBD \
                        --fastqc_args "--outdir /shared_lab/20180226_RNAseq_2017OAExp/DNAm/20190719_trim_galore_files --threads 18" \
                        $i \
                        $j \
                        2> /shared_lab/20180226_RNAseq_2017OAExp/DNAm/20190719_fastqc_trim_10bp_Cvirginica_MBD/stderr.log
                fi
        done
done
```
[File Link](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/DNAm/scripts/01_seq-quality-trim.sh)

**Outputs**
* Path: `/shared_lab/20180226_RNAseq_2017OAExp/DNAm/processed_samples/`
  * Folder `01_trimmed/20190719_fastqc_trim_10bp_Cvirginica_MBD`
      * Example File: `17162_DNAm_R1_val_1.fq.gz`

### Step Four - Preparing genome with Bismark for sample mapping

A bisulfite converted genome was prepared using the current *C. virginica* genome available NCBI using the `bismark_genome_preparation` command in `bismark`:  and `hisat2` as the aligner (it's splice aware compared to bowtie).

Code for bowtie2:
```
bismark_genome_preparation --bowtie2 --genomic_composition --parallel 10 --verbose /shared_lab/20180226_RNAseq_2017OAExp/DNAm/reference/genome > bismark_genomePrepartion_bowtie2_log.txt
```

Code for hisat2 (splice aware aligner):
```
bismark_genome_preparation --hisat2 --genomic_composition --path_to_aligner /home/downey-wall.a/software/hisat2-2.1.0 --parallel 10 --verbose /shared_lab/20180226_RNAseq_2017OAExp/DNAm/reference/genome > bismark_genomePrepartion_hisat2_log.txt
```

**Outputs**
* Path : `/shared_lab/20180226_RNAseq_2017OAExp/DNAm/reference/bsgenome/`
    * Hisat2 folder : `wHiSAT`
    * Bowtie2 folder : `wBowtie2`
      * `bismark_genomePrepartion_log.txt` : Command output log when using the `--verbose` flag.
      * `genomic_nucleotide_frequencies.txt` :  Text file with occurances of nucleotides in genome
      * `Bisulfite_Genome` : Folder with conversion information and locations

### Step Five - Aligning Reads
The trimmed reads were mapped to the bisulfite treated reference genome (created in the previous step) in order to determine the raw counts of methylated to unmethylated cytosines at each locus. This step was saved as several outputs including a sorted `.bam` file and a compressed `.txt` file with each row as a unique CpG.

Step one, mapping with bismark:
```
#!/bin/bash

trimmed="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/singleSample_trimScript_20190719/20190719_fastqc_trim_10bp_Cvirginica_MBD"
output="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/singleSample_trimScript_20190719/bismark/bowtie2"
# genome="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/reference/bsgenome/wBowtie2/"

cd $output

counter=1

echo "Start sample # (1-24)"
read lower
echo "End sample # (1-24)"
read upper

for i in $( ls $trimmed/*R1*gz); do
        file=$( echo $i | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 1)

        if [ "$counter" -ge "$lower" ]
        then

            	if [ "$counter" -le "$upper" ]
                then

                    	echo "Sample $counter"
                        echo "File Path : $i"
                        echo "File Name : $file"

                        # Run Bismark
                        bismark --non_directional -p 2\
                        --score_min L,0,-0.8 \
                        /shared_lab/20180226_RNAseq_2017OAExp/DNAm/reference/bsgenome/wBowtie2/ \
                        -1 $trimmed/${file}*R1*gz \
                        -2 $trimmed/${file}*R2*gz \
                        -o $output

                        # Deduplicate removal
                        deduplicate_bismark -p --bam \
                        $output/${file}*bam \
                        --output_dir $output

                        # Samtools reorder
                        samtools sort $output/${file}*deduplicated.bam \
                        -o $output/${file}dedup.sorted.bam

                        # Methylation Extraction
                        bismark_methylation_extractor -p --bedGraph --scaffolds --counts ${file}*deduplicated.bam --multicore 20

                        # Bismark2bedGraph
                        bismark2bedGraph --zero_based CpG_CTOB_${file}_DNAm_R1_val_1_bismark_bt2_pe.deduplicated.txt -o ${file}_CTOB
                        bismark2bedGraph --zero_based CpG_CTOT_${file}DNAm_R1_val_1_bismark_bt2_pe.deduplicated.txt -o ${file}_CTOT
                        bismark2bedGraph --zero_based CpG_OB_${file}_DNAm_R1_val_1_bismark_bt2_pe.deduplicated.txt -o ${file}_OB
                        bismark2bedGraph --zero_based CpG_OT_${file}_DNAm_R1_val_1_bismark_bt2_pe.deduplicated.txt -o ${file}_OT


                fi
        fi

	counter=$((counter+1))
done
```

**Outputs**
* Path: `/shared_lab/20180226_RNAseq_2017OAExp/DNAm/processed_samples/02_bismark`
   * Folder: `bowtie2`
       * `17145dedup.sorted.bam` : Sorted bam file with duplicates removed
       * `17145_DNAm_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz` : Text file with methy and ummethyl counts for each CpG (used for downstream analysis)

### Step Six - Cytosine Summary Reports
Summarizes the number of methylated and umethylated cytosine counts for each cytosine (in the context of the CpG motif).

Creating cytosine count files for each sample:
```
#!/bin/bash

dir="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/singleSample_trimScript_20190719/bismark/bowtie2"
genome="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/reference/genome"
out="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/singleSample_trimScript_20190719"

#cd $output
counter=1

echo "Start sample # (1-24)"
read lower
echo "End sample # (1-24)"
read upper

for i in $( ls $dir/*DNAm_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz); do
        file=$( echo $i | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 1)

        if [ "$counter" -ge "$lower" ]
        then

            	if [ "$counter" -le "$upper" ]
                then

                    	echo "Sample $counter"
                        echo "File Path : $i"
                        echo "File Name : $file"

                        coverage2cytosine ${i} \
                        --genome_folder $genome \
                        --dir $out \
                        -o ${file}_CytoSummary

                fi
        fi

	counter=$((counter+1))
done
```

**Output**
* Path: `/shared_lab/20180226_RNAseq_2017OAExp/DNAm/processed_samples`
	* Folder: `03_CytoSummaries`
		* `17090_CytoSummary.CpG_report.txt`

### Step Seven - Creating Count matrices (loci x sample)
Use custom r script to create a matrix for methylated and unmethylated cytosines, as well as a metaData sheet for storing information about each locus.

R script:
```
library(data.table)
library(dplyr)

sprintf("Reading in data and initializing methylation count matrix...")
# Set working directory to base folder where your raw files and bismark outputs reside
setwd("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/singleSample_trimScript_20190719")

# Specific folder where cov matrix are.
bis_outputs <- "CytoReports"
file_outputs <- "countMatrix"

files <- list.files(bis_outputs)
count_files <- files[grep("CytoSummary.CpG_report.txt",files)]
samples <- substr(count_files,1,5)

init <- fread(paste0(bis_outputs,"/",samples[1],"_CytoSummary.CpG_report.txt"),header=FALSE)

sum.table <- init[,c(1,2,3,6,7)]

init_C <- init[,4]
colnames(init_C) <- samples[1]
init_T <- init[,5]
colnames(init_T) <- samples[1]

sprintf("Starting counting....")
n<- 23
pb <- txtProgressBar(min = 0, max = n, style = 3)

for(i in 2:length(samples)){
  Sys.sleep(0.001)
  setTxtProgressBar(pb, i-1)

  temp <- fread(paste0(bis_outputs,"/",samples[i],"_CytoSummary.CpG_report.txt"),header=FALSE)
  temp_C <- temp[,4]
  temp_T <- temp[,5]

  colnames(temp_C) <- samples[i]
  colnames(temp_T) <- samples[i]

  init_C <- cbind(init_C,temp_C)
  init_T <- cbind(init_T,temp_T)
}
close(pb)

sprintf("Completed merging and reformatting matrices...")
#final_C <- data.frame(init_C[,-1])
#final_T <- data.frame(init_T[,-1])

#row.names(final_C) <- unlist(init_C[,1])
#row.names(final_T) <- unlist(init_T[,1])

#final_C[is.na(final_C)]<- 0
#final_T[is.na(final_T)]<- 0

#colnames(final_C) <- samples
#colnames(final_T) <- samples

final_S <- init_C + init_T
final_B<- init_C/final_S

sprintf("Saving files...")
fwrite(sum.table,paste0(file_outputs,"/All_CytoSum_summaryTable.csv"))
saveRDS(sum.table,paste0(file_outputs,"/All_CytoSum_summaryTable.RData"))


fwrite(init_C,paste0(file_outputs,"/All_CytoSum_methylCountMatrix.csv"))
fwrite(init_T,paste0(file_outputs,"/All_CytoSum_unmethylCountMatrix.csv"))
fwrite(final_S,paste0(file_outputs,"/All_CytoSum_TotalCountMatrix.csv"))
fwrite(final_B,paste0(file_outputs,"/All_CytoSum_BetaMatrix.csv"))

saveRDS(init_C,paste0(file_outputs,"/All_CytoSum_methylCountMatrix.RData"))
saveRDS(init_T,paste0(file_outputs,"/All_CytoSum_unmethylCountMatrix.RData"))
saveRDS(final_S,paste0(file_outputs,"/All_CytoSum_TotalCountMatrix.RData"))
saveRDS(final_B,paste0(file_outputs,"/All_CytoSum_BetaMatrix.RData"))
```

**Outputs**
* Path: `/shared_lab/20180226_RNAseq_2017OAExp/DNAm/processed_samples`
	* Folder: `04_countMatrix`
		* Example File: `All_CytoSum_BetaMatrix.RData`




