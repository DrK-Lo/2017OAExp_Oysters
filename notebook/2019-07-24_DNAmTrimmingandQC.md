# Performing trimming and QC on raw MBD-BSeq reads

Trimming and QC was performed on the concatenated read files using `trim_galore`. The parameters were selected based on the script Yaamini used for processing the 2016 oysters gonadal tissue experiment. This includes a 10 bp trim on either end of the read and the remaining parameters are left as default.

**Modifications to Original Script**
* I modified the original script created by Yaamini by using a `for` loop to cycle through each sample (R1 and R2) and perform the trim_galore command for each pair of reads.

**Potential Issues**
* I noticed that it is taking a long time to process each sample (it's been running for ~1.5-2 days and it hasn't finished processing all 24 samples). I am not sure if this normal, but seems like a long time for trimming and QC.
* Looking into the log files there appear to be some differences in the adapters being identified by `trim_galore`. This seems problematic given I would have expected them to be prepped all in the same way. 

**Script**
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

**Screen Shot of log outputs for two samples with different adapters identified**

