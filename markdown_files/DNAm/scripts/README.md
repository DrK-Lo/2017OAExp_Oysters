# Scripts


### Step One - Downloading the Data

Data was sent to the Roberts lab after sequencing and stored on their server [(Link for details)](). Below is the code used to download the data from their server:
```
FILL IN CODE 
```
Then moved to the appropriate working directory: `/shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp`

### Step Two - Concatenating the Data
Paired end sequencing for each sample was run across three lanes. These needed to be first concatenated into a single file and renamed so that the file ID matched the sample ID used elsewhere. This was done in a single step, using the `cat` function to combine files and storing the output in a new file that was renamed using the sampleID.

Raw file example: `zr2576_10_s1_R1.fastq.gz`
  * `zr2576_10` : Sample ID from the sequencer
  * `s1` : sequence lane (i.e. `s1`-`s3`)
  * `R1` : read 1 (paired end reads, `R1` or `R2`)
  
Bash script used to concatenate files:
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

Outputs were saved in the same directory. Example of file output: `17203_DNAm_R1.fastq.gz`

### Preparing genome with Bismark for sample mapping


