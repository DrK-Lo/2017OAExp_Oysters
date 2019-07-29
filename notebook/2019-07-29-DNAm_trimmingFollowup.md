# Follow up to potential DNAm trimming issues 

See [HERE](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/notebook/2019-07-24_DNAmTrimmingandQC.md) for initial discussion of observed issues using trim galore with coding loop.

### Follow up Overview
I reran the `trim_galore` script this time running it as a single command with all samples included. I ran it first without providing a specific flag for the adapter (auto detection mode) and also using the `--illumina` flag, since I expected the sequences to have been prepped with illumina adapters. I then compared the outputs to these two approached with those from the looped version of the script.

### Results
I found that they led to potentially different results, depending on the sample. When the single command was used and no adapter information was provided it treated each sample as if it had the nextera adapter. This is likely because the auto detection protocol determines the appropriate adapter based on the first 1 million reads. This means that it determined the adapter based on the first sample, `17005_DNAm_R1.fastq.gz`, processed. This is consistent with the single sample approach, which also found that the first sample had a nextera adapter. However, unlike the loop script, which tried to determine the appropriate adapter for each sample, the group command applied the nextera adapter to each sample. Using the `diff -s FILE1 FILE2` command, I saw that same output was produced regardless of script (single of loop) if they both identified the same adapter, but was different if they did not agree (i.e. nextera vs illumina).

### Current Outputs
I now have three versions of the trimming process
  1) Trimming performed for each sample with sample specific adapter detection (based on the loop script)
    * File Directory on Server : `/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/singleSample_trimScript_20190719`
  2) Trimming performed for each sample with single adapter detection based on first sample (single command script with auto detection)
    * File directory on server : `/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/groupSample_trimscipt_20190725`
  3) Trimming performed for each sample with forced adapter (illumina)
    * File directory on server (currently, will be moved to specific folder within the `trimmedSamples` folder when finished running) : `/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples`

### Next Steps
Need to confirm with Yaamini and Steven which of these three approaches seem like the correct way to go.
