Count Comparison
================
adowneywall
6/8/2019

### Data

``` r
STAR_unstrand <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/STAR_pipeline/rawCounts/rerun3_GeneCountMatrix_unstranded.rds")
STAR_strand <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/STAR_pipeline/rawCounts/rerun3_GeneCountMatrix_Stranded.rds")
STAR_reverse <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/STAR_pipeline/rawCounts/rerun3_GeneCountMatrix_Reverse.rds")
Salmon <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/salmon_pipeline/run20180512_gene_countMatrix_.RData")
GeneCounts <- read.delim("~/Github/2017OAExp_Oysters/input_files/RNA/C_virginica_gene_count_final.txt",header=TRUE,sep="",row.names=1)
```

### Count Summary

``` r
sal_match <- regmatches(row.names(Salmon),regexpr("LOC[0-9]*", row.names(Salmon)))
Salmon_red <- Salmon[match(sal_match,row.names(Salmon)),]

star_match <- regmatches(row.names(STAR_unstrand),regexpr("LOC[0-9]*", row.names(STAR_unstrand)))
STAR_und_red <- STAR_unstrand[match(star_match,row.names(STAR_unstrand)),]
STAR_str_red <- STAR_strand[match(star_match,row.names(STAR_strand)),]
STAR_rev_red <- STAR_reverse[match(star_match,row.names(STAR_reverse)),]
```

**Total counts for each approach**

``` r
# Counts based on Salmon Mapper
sum(Salmon_red)
```

    ## [1] 437419764

``` r
## Counts based on STAR --quantMode

# unstranded Counts (column 2)
sum(STAR_und_red)
```

    ## [1] 401635077

``` r
# 1st read Counts (column 3)
sum(STAR_str_red)
```

    ## [1] 3932439

``` r
# 2nd read Counts (column 4)
sum(STAR_reverse)
```

    ## [1] 399654362

``` r
# Gene counts using STAR original mapping and HT-Seq (-s flag which is equivalent to 1st read count in STAR --quantMode)
sum(GeneCounts)
```

    ## [1] 4742788

### Simple correlation between different mapping approaches and strand choice

``` r
# For the sake of comparisons I reduced all count matrices from each approach to only genes shared among all methods.
# Also for this simple test for water ever reason the final sample did not run with the STAR script so sample 24 was removed from comparison in both teh GeneCount and Salmon matrices
GeneCount_match <- GeneCounts[(match(row.names(Salmon_red),row.names(GeneCounts))),]
GeneCount_match <- GeneCount_match[,1:23]
STAR_unstr_match <- STAR_und_red[(match(row.names(Salmon_red),row.names(STAR_und_red))),]
STAR_rev_match <- STAR_rev_red[(match(row.names(Salmon_red),row.names(STAR_rev_red))),]
STAR_str_match <- STAR_str_red[(match(row.names(Salmon_red),row.names(STAR_str_red))),]
Salmon_red <- Salmon_red[,1:23]

#Correlations (mean and variance)
corCompare <- function(x,y){
  temp <- NULL
  for(i in 1:ncol(x)){
    temp <- c(temp,cor.test(x[,i],y[,i])$estimate)
  }
  out <- c(mean(temp),var(temp))
  return(out)
}

S1_GC <- corCompare(STAR_str_match,GeneCount_match)
S1_SU <- corCompare(STAR_str_match,STAR_unstr_match)
S1_S2 <- corCompare(STAR_str_match,STAR_rev_match)
S1_Sal <- corCompare(STAR_str_match,Salmon_red)

S2_GC <- corCompare(STAR_rev_match,GeneCount_match)
S2_SU <- corCompare(STAR_rev_match,STAR_unstr_match)
S2_Sal <- corCompare(STAR_rev_match,Salmon_red)

SU_GC <- corCompare(STAR_unstr_match,GeneCount_match)
SU_Sal <- corCompare(STAR_unstr_match,Salmon_red)

Sal_GC <- corCompare(Salmon_red,GeneCount_match)

corMat <- matrix(ncol=5,nrow=5)
nameMat <- c("STAR_Unstranded","STAR_Strand1","STAR_Strand2","STAR_HtSeq","Salmon")
row.names(corMat) <- nameMat
colnames(corMat) <- nameMat
corMat[1,] <- c(1,S1_SU[1],S2_SU[1],SU_GC[1],SU_Sal[1])
corMat[2,] <- c(S1_SU[2],1,S1_S2[1],S1_GC[1],S1_Sal[1])
corMat[3,] <- c(S2_SU[2],S1_S2[2],1,S2_GC[1],S2_Sal[1])
corMat[4,] <- c(SU_GC[2],S1_GC[2],S2_GC[2],1,Sal_GC[1])
corMat[5,] <- c(SU_Sal[2],S1_Sal[2],S1_Sal[2],Sal_GC[2],1)

library(kableExtra)
kable(corMat) %>% kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

STAR\_Unstranded

</th>

<th style="text-align:right;">

STAR\_Strand1

</th>

<th style="text-align:right;">

STAR\_Strand2

</th>

<th style="text-align:right;">

STAR\_HtSeq

</th>

<th style="text-align:right;">

Salmon

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

STAR\_Unstranded

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

0.1507156

</td>

<td style="text-align:right;">

0.9997707

</td>

<td style="text-align:right;">

0.0795862

</td>

<td style="text-align:right;">

0.5226464

</td>

</tr>

<tr>

<td style="text-align:left;">

STAR\_Strand1

</td>

<td style="text-align:right;">

0.0012083

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

0.1318888

</td>

<td style="text-align:right;">

0.9784128

</td>

<td style="text-align:right;">

0.0862911

</td>

</tr>

<tr>

<td style="text-align:left;">

STAR\_Strand2

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0013365

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

0.0609335

</td>

<td style="text-align:right;">

0.5226677

</td>

</tr>

<tr>

<td style="text-align:left;">

STAR\_HtSeq

</td>

<td style="text-align:right;">

0.0002316

</td>

<td style="text-align:right;">

0.0000877

</td>

<td style="text-align:right;">

0.0001945

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

0.1046505

</td>

</tr>

<tr>

<td style="text-align:left;">

Salmon

</td>

<td style="text-align:right;">

0.0062442

</td>

<td style="text-align:right;">

0.0002448

</td>

<td style="text-align:right;">

0.0002448

</td>

<td style="text-align:right;">

0.0004088

</td>

<td style="text-align:right;">

1.0000000

</td>

</tr>

</tbody>

</table>

### **Conclusions**

Here I compare the outputs from three different mapping approaches that
use two different RNA read mappers (STAR and Salmon). In addition, I
examine the difference between the counts produced when using the
`--quantMode` flag in STAR vs. using HT-Seq with the `-s` flag to
generate counts after the STAR mapping step. Note, I compared all count
outputs from the STAR –quantMode\` approach, these included;
**unstranded**, **first strand**, and **second strand** counts.

**Brief details about the different runs for each approach**

  - STAR with `--quanMode`: default settings were used to map reads.
    See:
    <https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/notebook/2019-06-04_STARrerunWithKMcodeComparison.md>
    for sample code.

  - STAR with `HT-Seq -s` counts: used slightly relaxed mapping quality
    (`--outSAMmapqUnique 40` vs 60), and non default criteria for
    `--outFilterMatchNminOverLread 0.17`and
    `--outFilterScoreMinOverLread 0.17`. See for full description:
    <https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/STAR_pipeline/outdated/brett_script.md>.

  - Salmon run with default values. See:
    <https://github.com/epigeneticstoocean/2017OAExp_Oysters/tree/master/markdown_files/Salmon_pipeline>
    for code.

**Major Thoughts**

1)  The unstranded and 2nd strand counts generated from the `--quanMode`
    flag are very similar. Generally the 2nd strand has consistently
    more reads than the 1st strand (sum is several order of magnitude
    larger), however this is not absolute. Further reading suggests that
    the unstranded counts may seem like a sum of the 1st and 2nd strand,
    but this is not necessary true (this can be see in the fact that
    unstranded and 1st strand are not particularly well correlated).

2)  The first strand counts are very well correlated with the `HT-Seq
    -s` counting approach. This makes sense because in theory the
    quantification of the first strand with the `quantMode` flag in STAR
    should be equivalent to the `-s` flag in HT-Seq. The reason they are
    not identical is likely more to do with the slightly different input
    parameters in the two different STAR runs used to generate the count
    matrices.

3)  The Salmon run correlates the best with the unstranded or 2nd strand
    matrices, this makes some sense since I don’t think it really
    considers “strandedness” in the same way that STAR does. As a
    result, I might expect the unstranded counts to be most analogous.
    Furthermore, the relatively marginally agreement between these
    methods may reflect some pretty fundamental differences between the
    two approaches. Not sure yet, what these means for the inferences
    between these two methods.

**Other outcomes**

This comparison gave me some important insight into why we were seeing
drastically different results from the original STAR run Brett performed
using HT-Seq fro the gene counts and what I got from Salmon. Despite
still only marginal agreement in the genes between STAR and Salmon, the
massive loss of counts had everything to do with the strand being
considered. Brett selected the `-s` flage in HT-Seq which only considers
the first read, which for reasons im still thinking about was
significantly smaller (across all genes) than the second strand. Also I
don’t see a reason to continue using Ht-Seq at the moment given that
`--quantMode` within STAR performs similarly and gives the flexibility
of considering one or all of the strand options post mapping.
