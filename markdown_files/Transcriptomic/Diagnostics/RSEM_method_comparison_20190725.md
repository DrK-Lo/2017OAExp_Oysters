RSEM Method Comparison
================
adowneywall
7/26/2019

## Overview

I was having trouble getting the RSEM transcript quantification command
to work properly with STAR produced bam files. The error created
mentioned that there were unequal transcripts identified in the STAR bam
file vs. the transcriptomic index created for RSEM (I used the star
mapper to create this index initially, and used the same GTF files to
produce both the STAR and RSEM index folders).

To fix this I decided to try three things.

  - 1)  Use a GTF file that only uses transcripts identified by GNOMON.
        This excludes ~ 24-50 identified by refSeq. These tended to be
        poorly annotated anyways.

  - 2)  I tried running STAR using either an index created using the
        STAR command for generating a genomic index or the RSEM command
        with the star flag.

  - 3)  I used the one step RSEM command with the RSEM index using the
        same GTF described above to perform the STAR mapping and RSEM
        quantification in a single command.

**Goal of this script is to compared the outputs of three RSEM
quantification strategies**

1)  bam\_STAR : STAR (w STAR Index) -\> RSEM
2)  bam\_RSEM : STAR (w RSEM Index) -\> RSEM
3)  rsem : RSEM (one step
option)

## Data

``` r
## Results from running RSEM quantification on STAR produced bam files using an STAR genome index for the STAR mapping
bam_star_genes <- read.delim("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/extra/17005_RSEM_test/RSEM_fromBAM_withSTARIndex.genes.results")
bam_star_isoforms <-  read.delim("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/extra/17005_RSEM_test/RSEM_fromBAM_withSTARIndex.isoforms.results")

## Results from running RSEM quantification on STAR produced bam files using an RSEM genome index for the STAR mapping
bam_rsem_genes <- read.delim("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/extra/17005_RSEM_test/RSEM_fromBAM_withRSEMIndex.genes.results")
bam_rsem_isoforms <- read.delim("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/extra/17005_RSEM_test/RSEM_fromBAM_withRSEMIndex.isoforms.results")

## Results from one step (STAR and RSEM in single command) run
rsem_genes <- read.delim("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/extra/17005_RSEM_test/RSEM_oneStep_withRSEMIndex.genes.results")
rsem_isoforms <- read.delim("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/extra/17005_RSEM_test/RSEM_oneStep_withRSEMIndex.isoforms.results")
```

### Pearsons Correlation

``` r
# Genes
corr_df <- data.frame(Output=colnames(bam_star_genes)[3:7],Gene_Feature="Gene",
                      bamStar_RSEM_corr=0,bamStar_RSEM_corr_pval=0,
                      bamStar_bamRSEM_corr=0,bamStar_bamRSEM_pval=0,
                      RSEM_bamRSEM_corr=0,RSEM_bamRSEM_pval=0)
for (i in 3:ncol(bam_star_genes)){
  corr_df[i-2,3] <- cor.test(bam_star_genes[,i],rsem_genes[,i])$estimate
  corr_df[i-2,4] <- cor.test(bam_star_genes[,i],rsem_genes[,i])$p.value
  corr_df[i-2,5] <- cor.test(bam_star_genes[,i],bam_rsem_genes[,i])$estimate
  corr_df[i-2,6] <- cor.test(bam_star_genes[,i],bam_rsem_genes[,i])$p.value
  corr_df[i-2,7] <- cor.test(rsem_genes[,i],bam_rsem_genes[,i])$estimate
  corr_df[i-2,8] <- cor.test(rsem_genes[,i],bam_rsem_genes[,i])$p.value
}
# Isoforms
corr_df_iso <- data.frame(Output=colnames(bam_star_isoforms)[3:7],Gene_Feature="Isoform",
                      bamStar_RSEM_corr=0,bamStar_RSEM_corr_pval=0,
                      bamStar_bamRSEM_corr=0,bamStar_bamRSEM_pval=0,
                      RSEM_bamRSEM_corr=0,RSEM_bamRSEM_pval=0)

for (i in 4:ncol(bam_star_isoforms)){
  corr_df_iso[i-3,3] <- cor.test(bam_star_isoforms[,i],rsem_isoforms[,i])$estimate
  corr_df_iso[i-3,4] <- cor.test(bam_star_isoforms[,i],rsem_isoforms[,i])$p.value
  corr_df_iso[i-3,5] <- cor.test(bam_star_isoforms[,i],bam_rsem_isoforms[,i])$estimate
  corr_df_iso[i-3,6] <- cor.test(bam_star_isoforms[,i],bam_rsem_isoforms[,i])$p.value
  corr_df_iso[i-3,7] <- cor.test(rsem_isoforms[,i],bam_rsem_isoforms[,i])$estimate
  corr_df_iso[i-3,8] <- cor.test(rsem_isoforms[,i],bam_rsem_isoforms[,i])$p.value
}

# Merging dataframes
corr_df_whole <- rbind.data.frame(corr_df,corr_df_iso)

kableExtra::kable(corr_df_whole) %>% kableExtra::kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Output

</th>

<th style="text-align:left;">

Gene\_Feature

</th>

<th style="text-align:right;">

bamStar\_RSEM\_corr

</th>

<th style="text-align:right;">

bamStar\_RSEM\_corr\_pval

</th>

<th style="text-align:right;">

bamStar\_bamRSEM\_corr

</th>

<th style="text-align:right;">

bamStar\_bamRSEM\_pval

</th>

<th style="text-align:right;">

RSEM\_bamRSEM\_corr

</th>

<th style="text-align:right;">

RSEM\_bamRSEM\_pval

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

length

</td>

<td style="text-align:left;">

Gene

</td>

<td style="text-align:right;">

0.9997856

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.9997856

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

effective\_length

</td>

<td style="text-align:left;">

Gene

</td>

<td style="text-align:right;">

0.9997855

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.9997855

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

expected\_count

</td>

<td style="text-align:left;">

Gene

</td>

<td style="text-align:right;">

0.9995225

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.9995225

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

TPM

</td>

<td style="text-align:left;">

Gene

</td>

<td style="text-align:right;">

0.9988203

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.9988203

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

FPKM

</td>

<td style="text-align:left;">

Gene

</td>

<td style="text-align:right;">

0.9988203

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.9988203

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

length

</td>

<td style="text-align:left;">

Isoform

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

effective\_length

</td>

<td style="text-align:left;">

Isoform

</td>

<td style="text-align:right;">

0.9993526

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.9993526

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

expected\_count

</td>

<td style="text-align:left;">

Isoform

</td>

<td style="text-align:right;">

0.9986665

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.9986665

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

TPM

</td>

<td style="text-align:left;">

Isoform

</td>

<td style="text-align:right;">

0.9986665

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.9986665

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:left;">

FPKM

</td>

<td style="text-align:left;">

Isoform

</td>

<td style="text-align:right;">

0.9755454

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.9755454

</td>

<td style="text-align:right;">

0

</td>

</tr>

</tbody>

</table>
