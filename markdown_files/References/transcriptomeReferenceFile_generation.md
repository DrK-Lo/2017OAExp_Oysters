---
title: "Script for creating reference transcriptomes"
author: "adowneywall"
date: "5/12/2019"
output:
  html_document:
     keep_md: true
editor_options: 
  chunk_output_type: console
---



**Description**: This script takes the latest version of the the Eastern oyster genome and splits in into a data.frame.  

**Reference files**  
[Click here for NCBI ftp URL](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/GCF_002022765.2_C_virginica-3.0)  

**Original File IDs** 
* ```GCF_002022765.2_C_virginica-3.0_rna.fna```  
* ```GCF_002022765.2_C_virginica-3.0_rna_from_genomic.fna```  

**Saved dataframe files**  
* ```transcriptome_table.RData```  
* ```transcriptome_fromGenome_table.RData```  

### Reading in a transcriptome data  

```r
trans <- readLines("/home/downeyam/Github/salmon_tutorial/GCF_002022765.2_C_virginica-3.0_rna.fna")
save_trans <- grep(">",trans)
trans_select <- trans[save_trans]

## Split up trans file
# Full ID
t_fullID <- sub(">(.+?)\\s.*", "\\1", trans_select)
#Gene Prediction
t_Predict <- sub(".*PREDICTED:\\s(.*)\\s\\(.*", "\\1", trans_select)
# Gene location (matches with other file)
t_location <- sub(".*\\((.*)\\).*","\\1",trans_select)
# Gene Variant
t_transVariant <- sub(".*,\\s(.*),.*","\\1",trans_select)
# Gene Type
t_type <- sub(".*,\\s(.*)","\\1",trans_select)

t_rna <- data.frame(fullID=t_fullID,location=t_location,transVariant=t_transVariant,predict=t_Predict,type=t_type)
#saveRDS(t_rna,"/home/downeyam/Github/2017OAExp_Oysters/input_files/transcriptome_table.RData")

kable(head(t_rna)) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> fullID </th>
   <th style="text-align:left;"> location </th>
   <th style="text-align:left;"> transVariant </th>
   <th style="text-align:left;"> predict </th>
   <th style="text-align:left;"> type </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> XM_022430339.1 </td>
   <td style="text-align:left;"> LOC111138521 </td>
   <td style="text-align:left;"> transcript variant X3 </td>
   <td style="text-align:left;"> Crassostrea virginica purine nucleoside phosphorylase-like </td>
   <td style="text-align:left;"> mRNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> XM_022430340.1 </td>
   <td style="text-align:left;"> LOC111099031 </td>
   <td style="text-align:left;"> transcript variant X1 </td>
   <td style="text-align:left;"> Crassostrea virginica uncharacterized LOC111099031 </td>
   <td style="text-align:left;"> mRNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> XM_022430341.1 </td>
   <td style="text-align:left;"> LOC111099031 </td>
   <td style="text-align:left;"> transcript variant X2 </td>
   <td style="text-align:left;"> Crassostrea virginica uncharacterized LOC111099031 </td>
   <td style="text-align:left;"> mRNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> XM_022430342.1 </td>
   <td style="text-align:left;"> LOC111099031 </td>
   <td style="text-align:left;"> transcript variant X3 </td>
   <td style="text-align:left;"> Crassostrea virginica uncharacterized LOC111099031 </td>
   <td style="text-align:left;"> mRNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> XM_022430343.1 </td>
   <td style="text-align:left;"> LOC111099032 </td>
   <td style="text-align:left;"> &gt;XM_022430343.1 PREDICTED: Crassostrea virginica 5-hydroxytryptamine receptor-like (LOC111099032), mRNA </td>
   <td style="text-align:left;"> Crassostrea virginica 5-hydroxytryptamine receptor-like </td>
   <td style="text-align:left;"> mRNA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> XM_022430344.1 </td>
   <td style="text-align:left;"> LOC111099034 </td>
   <td style="text-align:left;"> &gt;XM_022430344.1 PREDICTED: Crassostrea virginica uncharacterized LOC111099034 (LOC111099034), mRNA </td>
   <td style="text-align:left;"> Crassostrea virginica uncharacterized LOC111099034 </td>
   <td style="text-align:left;"> mRNA </td>
  </tr>
</tbody>
</table>

### Transcriptomic data from genome  

```r
trans_fg <- readLines("/home/downeyam/Github/salmon_tutorial/GCF_002022765.2_C_virginica-3.0_rna_from_genomic.fna")
save_trans_fg <- grep(">",trans_fg)
trans_fg_select <- trans_fg[save_trans_fg]

## Split up trans_fg file
#Gene Location
fg_location <- sub(".*gene=(LOC[0-9]+).*", "\\1", trans_fg_select)
# Gene ID
fg_geneID <- sub(".*GeneID:(.+?)].*", "\\1", trans_fg_select)
# Product
fg_product <- sub(".*product=(.+?)].*", "\\1", trans_fg_select)
# Transcript ID
fg_transcriptID <- sub(".*transcript_id=(.+?)].*", "\\1", trans_fg_select)
# Full ID
fg_fullID <- sub(">(.+?)\\s.*", "\\1", trans_fg_select)

fg <- data.frame(fullID=fg_fullID,location=fg_location,geneID=fg_geneID,product=fg_product,transcriptID=fg_transcriptID)
#saveRDS(fg,"/home/downeyam/Github/2017OAExp_Oysters/input_files/transcriptome_fromGenome_table.RData")
kable(head(fg)) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> fullID </th>
   <th style="text-align:left;"> location </th>
   <th style="text-align:left;"> geneID </th>
   <th style="text-align:left;"> product </th>
   <th style="text-align:left;"> transcriptID </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> lcl|NC_035780.1_ncrna_XR_002636969.1_1 </td>
   <td style="text-align:left;"> LOC111116054 </td>
   <td style="text-align:left;"> 111116054 </td>
   <td style="text-align:left;"> uncharacterized LOC111116054 </td>
   <td style="text-align:left;"> XR_002636969.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lcl|NC_035780.1_mrna_XM_022471938.1_2 </td>
   <td style="text-align:left;"> LOC111126949 </td>
   <td style="text-align:left;"> 111126949 </td>
   <td style="text-align:left;"> UNC5C-like protein </td>
   <td style="text-align:left;"> XM_022471938.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lcl|NC_035780.1_mrna_XM_022447324.1_3 </td>
   <td style="text-align:left;"> LOC111110729 </td>
   <td style="text-align:left;"> 111110729 </td>
   <td style="text-align:left;"> FMRFamide receptor-like, transcript variant X1 </td>
   <td style="text-align:left;"> XM_022447324.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lcl|NC_035780.1_mrna_XM_022447333.1_4 </td>
   <td style="text-align:left;"> LOC111110729 </td>
   <td style="text-align:left;"> 111110729 </td>
   <td style="text-align:left;"> FMRFamide receptor-like, transcript variant X2 </td>
   <td style="text-align:left;"> XM_022447333.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lcl|NC_035780.1_mrna_XM_022449924.1_5 </td>
   <td style="text-align:left;"> LOC111112434 </td>
   <td style="text-align:left;"> 111112434 </td>
   <td style="text-align:left;"> homeobox protein Hox-B7-like </td>
   <td style="text-align:left;"> XM_022449924.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lcl|NC_035780.1_mrna_XM_022461698.1_6 </td>
   <td style="text-align:left;"> LOC111120752 </td>
   <td style="text-align:left;"> 111120752 </td>
   <td style="text-align:left;"> ribulose-phosphate 3-epimerase-like </td>
   <td style="text-align:left;"> XM_022461698.1 </td>
  </tr>
</tbody>
</table>
