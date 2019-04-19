# Downloaded newest GFF files from NCBI and performed comparison with old files using gffcompare

### **New Programs Installed on cluster**  

`gffread` and `gffcompare` - http://ccb.jhu.edu/software/stringtie/gff.shtml

* `gffread` is can be used to validate, filter, convert and perform various other operations on GFF files (use gffread -h to see the various usage options). 
* `gffcompare` can be used to compare, merge, annotate and estimate accuracy of one or more GFF files (the "query" files), when compared with a reference annotation (also provided as GFF/GTF).

**Note**: I created an alias for the primary executables for both of these programs on my local ~/.bashrc on the comp5 cluster. They can be called by using `gffread_run` or `gffcompare_run`, respectively.


### Creating new GTF file from latest .gff file from NCBI
Oyster genome page : https://www.ncbi.nlm.nih.gov/genome/?term=eastern+oyster
Lastest .gff : ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/GCF_002022765.2_C_virginica-3.0/GCF_002022765.2_C_virginica-3.0_genomic.gff.gz

Converted .gff into .gtf using `gffread` in this line of code:  
` gffread_run GCF_002022765.2_C_virginica-3.0_genomic.gff -T -o top_v2.gtf`

The new file is within the ~/Genome/GFF directory and a link (`ln`) was made to a file of the same name in the ~/RNA directory.

### Next I compare the old .gtf with the new one

This line of code:  
`gffcompare_run -R -r top_v2.gtf -o refGTF top.gtf`

*Note you need to be in the RNA directory to do this*

This created several output summaries that compared the degree of overlap between the new (`top_v2.gtf`) and old (`top.gtf') .gtf files

**Summary Stats**  
```
Summary for dataset: top.gtf
Query mRNAs :   67891 in   39152 loci  (65422 multi-exon transcripts)
            (11436 multi-transcript loci, ~1.7 transcripts per locus)
Reference mRNAs :   67853 in   39152 loci  (65384 multi-exon)
Super-loci w/ reference transcripts:    39152
 -----------------| Sensitivity | Precision  |
        Base level:   100.0     |   100.0    |
        Exon level:   100.0     |   100.0    |
      Intron level:   100.0     |   100.0    |
Intron chain level:   100.0     |    99.9    |
  Transcript level:   100.0     |    99.9    |
       Locus level:   100.0     |   100.0    |

     Matching intron chains:   65384
       Matching transcripts:   67853
              Matching loci:   39152

          Missed exons:       0/352731  (  0.0%)
           Novel exons:       0/352757  (  0.0%)
        Missed introns:       0/310704  (  0.0%)
         Novel introns:       0/310704  (  0.0%)
           Missed loci:       0/39152   (  0.0%)
            Novel loci:       0/39152   (  0.0%)
Total union super-loci across all input datasets: 39152 
67891 out of 67891 consensus transcripts written in refGTF.annotated.gtf (0 discarded as redundant)
```

Full results can be found in a new folder called `gffcompare_out` in the RNA directory.
  
**Conclusions** :  
* There doesn't seem to be much of a significant difference between the two files, indicating the old .gtf file is still good. 
* I will stick with the old file for future use because Brett went in a manual added transcript IDs to this file for the rare cases that an entry was missing an id.
* Suggests another reason for small gene count matrix.
