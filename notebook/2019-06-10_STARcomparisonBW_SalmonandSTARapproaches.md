# Salmon vs. Star Mapping approaches

**Overview**: I looked at the mapping outcomes when you use the `--quantMode` option in STAR. Specifically, I compared the different count ouputs this option creates against the counts produced from our original STAR run (using HT-Seq to count expression) and the outputs from Salmon.

**Outcome**: Major outcomes include a better understanding of why we have substantially fewer (orders of magnitude fewer) counts when we performed the STAR w/ Ht_Seq option (it has to do with which strand is being included in the count). In addition, I showed that `STAR -quantMode` (with the correct strand) has high agreement with the HT-Seq based gene quantification and Salmon shows higher correlation with unstranded gene counts from STAR, but due to differences in their approach still only exhibit a correlation of ~0.5.

See Here for More details:
https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/extra/STAR_Salmon_Counting_Comparison.md
