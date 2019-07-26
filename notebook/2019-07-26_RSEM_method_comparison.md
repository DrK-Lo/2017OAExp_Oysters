## Compared Transcript Quantification from multiple STAR/RSEM mapping and quantification approaches

## Overview
I was having trouble getting the RSEM transcript quantification command to work properly with STAR produced bam files. The error created mentioned that there were unequal transcripts identified in the STAR bam file vs. the transcriptomic index created for RSEM (I used the star mapper to create this index initially, and used the same GTF files to produce both the STAR and RSEM index folders). 

To fix this I decided to try three things. 

* 1) Use a GTF file that only uses transcripts identified by GNOMON. This excludes ~ 24-50 identified by refSeq. These tended to be poorly annotated anyways.

* 2) I tried running STAR using either an index created using the STAR command for generating a genomic index or the RSEM command with the star flag. 

* 3) I used the one step RSEM command with the RSEM index using the same GTF described above to perform the STAR mapping and RSEM quantification in a single command.

![]()
