---
title: "STAR_Salmon_Mapping_Comparison"
author: "adowneywall"
date: "5/14/2019"
output: 
  html_document: 
    keep_md: true
editor_options: 
  chunk_output_type: console
---



### **Data**

```r
sa_c <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/salmon_RNA/run20180512_gene_countMatrix_.RData")  
sa_c <- sa_c[-c(1:205),] #removing some no LOC genes (most trna)  

st_c <- read.delim("/home/downeyam/Github/2017OAExp_Oysters/results/C_virginica_gene_count_final.txt",sep = " ")
```

### **Looking at non-overlapping transcripts**

```r
# Genes that are unique to the Salmon mapping
length(setdiff(row.names(sa_c),row.names(st_c)))
```

```
## [1] 0
```

```r
# None

# Genes that are unique to the STAR mapping
length(setdiff(row.names(st_c),row.names(sa_c)))
```

```
## [1] 836
```

```r
# 836 Genes

## Create dataframe of ordered genes unique to STAR mapping
st_genes <- setdiff(row.names(st_c),row.names(sa_c))
st_reduce <- st_c[match(st_genes,row.names(st_c)),]
sum_vec <- rowSums(st_reduce)
st_reduce <- st_reduce[rev(order(sum_vec)),]
```

Summary Table of Non-overlapping Transcripts Order by Count

```r
kable(st_reduce[1:20,]) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> RNA17005 </th>
   <th style="text-align:right;"> RNA17007 </th>
   <th style="text-align:right;"> RNA17013 </th>
   <th style="text-align:right;"> RNA17019 </th>
   <th style="text-align:right;"> RNA17069 </th>
   <th style="text-align:right;"> RNA17070 </th>
   <th style="text-align:right;"> RNA17072 </th>
   <th style="text-align:right;"> RNA17079 </th>
   <th style="text-align:right;"> RNA17090 </th>
   <th style="text-align:right;"> RNA17094 </th>
   <th style="text-align:right;"> RNA17099 </th>
   <th style="text-align:right;"> RNA17108 </th>
   <th style="text-align:right;"> RNA17122 </th>
   <th style="text-align:right;"> RNA17130 </th>
   <th style="text-align:right;"> RNA17142 </th>
   <th style="text-align:right;"> RNA17145 </th>
   <th style="text-align:right;"> RNA17162 </th>
   <th style="text-align:right;"> RNA17174 </th>
   <th style="text-align:right;"> RNA17176 </th>
   <th style="text-align:right;"> RNA17178 </th>
   <th style="text-align:right;"> RNA17181 </th>
   <th style="text-align:right;"> RNA17203 </th>
   <th style="text-align:right;"> RNA17211 </th>
   <th style="text-align:right;"> RNA17213 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> unknown_gene2 </td>
   <td style="text-align:right;"> 1079 </td>
   <td style="text-align:right;"> 820 </td>
   <td style="text-align:right;"> 1923 </td>
   <td style="text-align:right;"> 817 </td>
   <td style="text-align:right;"> 2294 </td>
   <td style="text-align:right;"> 1232 </td>
   <td style="text-align:right;"> 1661 </td>
   <td style="text-align:right;"> 1363 </td>
   <td style="text-align:right;"> 1547 </td>
   <td style="text-align:right;"> 1168 </td>
   <td style="text-align:right;"> 1229 </td>
   <td style="text-align:right;"> 1515 </td>
   <td style="text-align:right;"> 2694 </td>
   <td style="text-align:right;"> 847 </td>
   <td style="text-align:right;"> 1197 </td>
   <td style="text-align:right;"> 1479 </td>
   <td style="text-align:right;"> 977 </td>
   <td style="text-align:right;"> 1634 </td>
   <td style="text-align:right;"> 969 </td>
   <td style="text-align:right;"> 839 </td>
   <td style="text-align:right;"> 774 </td>
   <td style="text-align:right;"> 2291 </td>
   <td style="text-align:right;"> 1360 </td>
   <td style="text-align:right;"> 2197 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LOC111123407 </td>
   <td style="text-align:right;"> 1045 </td>
   <td style="text-align:right;"> 876 </td>
   <td style="text-align:right;"> 537 </td>
   <td style="text-align:right;"> 1010 </td>
   <td style="text-align:right;"> 789 </td>
   <td style="text-align:right;"> 756 </td>
   <td style="text-align:right;"> 403 </td>
   <td style="text-align:right;"> 421 </td>
   <td style="text-align:right;"> 492 </td>
   <td style="text-align:right;"> 672 </td>
   <td style="text-align:right;"> 645 </td>
   <td style="text-align:right;"> 489 </td>
   <td style="text-align:right;"> 563 </td>
   <td style="text-align:right;"> 444 </td>
   <td style="text-align:right;"> 641 </td>
   <td style="text-align:right;"> 618 </td>
   <td style="text-align:right;"> 734 </td>
   <td style="text-align:right;"> 745 </td>
   <td style="text-align:right;"> 844 </td>
   <td style="text-align:right;"> 620 </td>
   <td style="text-align:right;"> 964 </td>
   <td style="text-align:right;"> 629 </td>
   <td style="text-align:right;"> 674 </td>
   <td style="text-align:right;"> 823 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> unknown_gene1 </td>
   <td style="text-align:right;"> 281 </td>
   <td style="text-align:right;"> 275 </td>
   <td style="text-align:right;"> 524 </td>
   <td style="text-align:right;"> 320 </td>
   <td style="text-align:right;"> 756 </td>
   <td style="text-align:right;"> 379 </td>
   <td style="text-align:right;"> 478 </td>
   <td style="text-align:right;"> 460 </td>
   <td style="text-align:right;"> 215 </td>
   <td style="text-align:right;"> 385 </td>
   <td style="text-align:right;"> 297 </td>
   <td style="text-align:right;"> 259 </td>
   <td style="text-align:right;"> 1600 </td>
   <td style="text-align:right;"> 246 </td>
   <td style="text-align:right;"> 313 </td>
   <td style="text-align:right;"> 316 </td>
   <td style="text-align:right;"> 203 </td>
   <td style="text-align:right;"> 611 </td>
   <td style="text-align:right;"> 283 </td>
   <td style="text-align:right;"> 203 </td>
   <td style="text-align:right;"> 214 </td>
   <td style="text-align:right;"> 622 </td>
   <td style="text-align:right;"> 408 </td>
   <td style="text-align:right;"> 1233 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LOC111128891 </td>
   <td style="text-align:right;"> 175 </td>
   <td style="text-align:right;"> 99 </td>
   <td style="text-align:right;"> 141 </td>
   <td style="text-align:right;"> 151 </td>
   <td style="text-align:right;"> 206 </td>
   <td style="text-align:right;"> 154 </td>
   <td style="text-align:right;"> 103 </td>
   <td style="text-align:right;"> 92 </td>
   <td style="text-align:right;"> 98 </td>
   <td style="text-align:right;"> 141 </td>
   <td style="text-align:right;"> 136 </td>
   <td style="text-align:right;"> 102 </td>
   <td style="text-align:right;"> 115 </td>
   <td style="text-align:right;"> 112 </td>
   <td style="text-align:right;"> 105 </td>
   <td style="text-align:right;"> 135 </td>
   <td style="text-align:right;"> 94 </td>
   <td style="text-align:right;"> 110 </td>
   <td style="text-align:right;"> 125 </td>
   <td style="text-align:right;"> 114 </td>
   <td style="text-align:right;"> 98 </td>
   <td style="text-align:right;"> 155 </td>
   <td style="text-align:right;"> 131 </td>
   <td style="text-align:right;"> 147 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trnaa-ugc </td>
   <td style="text-align:right;"> 83 </td>
   <td style="text-align:right;"> 117 </td>
   <td style="text-align:right;"> 123 </td>
   <td style="text-align:right;"> 122 </td>
   <td style="text-align:right;"> 116 </td>
   <td style="text-align:right;"> 126 </td>
   <td style="text-align:right;"> 154 </td>
   <td style="text-align:right;"> 73 </td>
   <td style="text-align:right;"> 166 </td>
   <td style="text-align:right;"> 90 </td>
   <td style="text-align:right;"> 125 </td>
   <td style="text-align:right;"> 124 </td>
   <td style="text-align:right;"> 121 </td>
   <td style="text-align:right;"> 79 </td>
   <td style="text-align:right;"> 78 </td>
   <td style="text-align:right;"> 130 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 218 </td>
   <td style="text-align:right;"> 180 </td>
   <td style="text-align:right;"> 106 </td>
   <td style="text-align:right;"> 76 </td>
   <td style="text-align:right;"> 68 </td>
   <td style="text-align:right;"> 152 </td>
   <td style="text-align:right;"> 151 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LOC111123439 </td>
   <td style="text-align:right;"> 111 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 148 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 272 </td>
   <td style="text-align:right;"> 156 </td>
   <td style="text-align:right;"> 149 </td>
   <td style="text-align:right;"> 154 </td>
   <td style="text-align:right;"> 156 </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 28 </td>
   <td style="text-align:right;"> 229 </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 116 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 97 </td>
   <td style="text-align:right;"> 207 </td>
   <td style="text-align:right;"> 58 </td>
   <td style="text-align:right;"> 29 </td>
   <td style="text-align:right;"> 142 </td>
   <td style="text-align:right;"> 103 </td>
   <td style="text-align:right;"> 127 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LOC111119537 </td>
   <td style="text-align:right;"> 123 </td>
   <td style="text-align:right;"> 123 </td>
   <td style="text-align:right;"> 116 </td>
   <td style="text-align:right;"> 42 </td>
   <td style="text-align:right;"> 49 </td>
   <td style="text-align:right;"> 79 </td>
   <td style="text-align:right;"> 89 </td>
   <td style="text-align:right;"> 87 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 128 </td>
   <td style="text-align:right;"> 112 </td>
   <td style="text-align:right;"> 70 </td>
   <td style="text-align:right;"> 93 </td>
   <td style="text-align:right;"> 69 </td>
   <td style="text-align:right;"> 48 </td>
   <td style="text-align:right;"> 113 </td>
   <td style="text-align:right;"> 112 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 124 </td>
   <td style="text-align:right;"> 55 </td>
   <td style="text-align:right;"> 43 </td>
   <td style="text-align:right;"> 46 </td>
   <td style="text-align:right;"> 120 </td>
   <td style="text-align:right;"> 87 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LOC111125281 </td>
   <td style="text-align:right;"> 53 </td>
   <td style="text-align:right;"> 84 </td>
   <td style="text-align:right;"> 96 </td>
   <td style="text-align:right;"> 57 </td>
   <td style="text-align:right;"> 124 </td>
   <td style="text-align:right;"> 83 </td>
   <td style="text-align:right;"> 65 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 91 </td>
   <td style="text-align:right;"> 64 </td>
   <td style="text-align:right;"> 89 </td>
   <td style="text-align:right;"> 67 </td>
   <td style="text-align:right;"> 92 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 81 </td>
   <td style="text-align:right;"> 54 </td>
   <td style="text-align:right;"> 91 </td>
   <td style="text-align:right;"> 112 </td>
   <td style="text-align:right;"> 61 </td>
   <td style="text-align:right;"> 67 </td>
   <td style="text-align:right;"> 49 </td>
   <td style="text-align:right;"> 80 </td>
   <td style="text-align:right;"> 61 </td>
   <td style="text-align:right;"> 67 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LOC111109071 </td>
   <td style="text-align:right;"> 59 </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 44 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 42 </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 83 </td>
   <td style="text-align:right;"> 87 </td>
   <td style="text-align:right;"> 85 </td>
   <td style="text-align:right;"> 163 </td>
   <td style="text-align:right;"> 172 </td>
   <td style="text-align:right;"> 78 </td>
   <td style="text-align:right;"> 72 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 101 </td>
   <td style="text-align:right;"> 110 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 62 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 64 </td>
   <td style="text-align:right;"> 62 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LOC111123854 </td>
   <td style="text-align:right;"> 60 </td>
   <td style="text-align:right;"> 138 </td>
   <td style="text-align:right;"> 94 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 131 </td>
   <td style="text-align:right;"> 78 </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:right;"> 63 </td>
   <td style="text-align:right;"> 78 </td>
   <td style="text-align:right;"> 71 </td>
   <td style="text-align:right;"> 86 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 136 </td>
   <td style="text-align:right;"> 78 </td>
   <td style="text-align:right;"> 57 </td>
   <td style="text-align:right;"> 99 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 85 </td>
   <td style="text-align:right;"> 47 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LOC111123567 </td>
   <td style="text-align:right;"> 277 </td>
   <td style="text-align:right;"> 48 </td>
   <td style="text-align:right;"> 65 </td>
   <td style="text-align:right;"> 67 </td>
   <td style="text-align:right;"> 38 </td>
   <td style="text-align:right;"> 69 </td>
   <td style="text-align:right;"> 47 </td>
   <td style="text-align:right;"> 59 </td>
   <td style="text-align:right;"> 53 </td>
   <td style="text-align:right;"> 60 </td>
   <td style="text-align:right;"> 43 </td>
   <td style="text-align:right;"> 76 </td>
   <td style="text-align:right;"> 60 </td>
   <td style="text-align:right;"> 61 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:right;"> 64 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 66 </td>
   <td style="text-align:right;"> 60 </td>
   <td style="text-align:right;"> 42 </td>
   <td style="text-align:right;"> 43 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trnat-ggu </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 74 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 97 </td>
   <td style="text-align:right;"> 113 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 129 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 85 </td>
   <td style="text-align:right;"> 52 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 65 </td>
   <td style="text-align:right;"> 82 </td>
   <td style="text-align:right;"> 63 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 69 </td>
   <td style="text-align:right;"> 48 </td>
   <td style="text-align:right;"> 145 </td>
   <td style="text-align:right;"> 149 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LOC111132625 </td>
   <td style="text-align:right;"> 110 </td>
   <td style="text-align:right;"> 63 </td>
   <td style="text-align:right;"> 41 </td>
   <td style="text-align:right;"> 49 </td>
   <td style="text-align:right;"> 46 </td>
   <td style="text-align:right;"> 61 </td>
   <td style="text-align:right;"> 71 </td>
   <td style="text-align:right;"> 92 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 49 </td>
   <td style="text-align:right;"> 86 </td>
   <td style="text-align:right;"> 43 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 61 </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 72 </td>
   <td style="text-align:right;"> 42 </td>
   <td style="text-align:right;"> 48 </td>
   <td style="text-align:right;"> 97 </td>
   <td style="text-align:right;"> 69 </td>
   <td style="text-align:right;"> 65 </td>
   <td style="text-align:right;"> 102 </td>
   <td style="text-align:right;"> 29 </td>
   <td style="text-align:right;"> 53 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LOC111105197 </td>
   <td style="text-align:right;"> 106 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 46 </td>
   <td style="text-align:right;"> 94 </td>
   <td style="text-align:right;"> 72 </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 55 </td>
   <td style="text-align:right;"> 65 </td>
   <td style="text-align:right;"> 49 </td>
   <td style="text-align:right;"> 50 </td>
   <td style="text-align:right;"> 55 </td>
   <td style="text-align:right;"> 47 </td>
   <td style="text-align:right;"> 46 </td>
   <td style="text-align:right;"> 46 </td>
   <td style="text-align:right;"> 73 </td>
   <td style="text-align:right;"> 67 </td>
   <td style="text-align:right;"> 49 </td>
   <td style="text-align:right;"> 45 </td>
   <td style="text-align:right;"> 61 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 73 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LOC111131090 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 70 </td>
   <td style="text-align:right;"> 76 </td>
   <td style="text-align:right;"> 125 </td>
   <td style="text-align:right;"> 45 </td>
   <td style="text-align:right;"> 70 </td>
   <td style="text-align:right;"> 28 </td>
   <td style="text-align:right;"> 67 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 89 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 46 </td>
   <td style="text-align:right;"> 57 </td>
   <td style="text-align:right;"> 78 </td>
   <td style="text-align:right;"> 58 </td>
   <td style="text-align:right;"> 88 </td>
   <td style="text-align:right;"> 68 </td>
   <td style="text-align:right;"> 53 </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 50 </td>
   <td style="text-align:right;"> 47 </td>
   <td style="text-align:right;"> 43 </td>
   <td style="text-align:right;"> 25 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LOC111122630 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 44 </td>
   <td style="text-align:right;"> 44 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 87 </td>
   <td style="text-align:right;"> 79 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 95 </td>
   <td style="text-align:right;"> 60 </td>
   <td style="text-align:right;"> 80 </td>
   <td style="text-align:right;"> 80 </td>
   <td style="text-align:right;"> 57 </td>
   <td style="text-align:right;"> 71 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 58 </td>
   <td style="text-align:right;"> 50 </td>
   <td style="text-align:right;"> 41 </td>
   <td style="text-align:right;"> 43 </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:right;"> 45 </td>
   <td style="text-align:right;"> 28 </td>
   <td style="text-align:right;"> 57 </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 73 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LOC111105411 </td>
   <td style="text-align:right;"> 43 </td>
   <td style="text-align:right;"> 88 </td>
   <td style="text-align:right;"> 54 </td>
   <td style="text-align:right;"> 123 </td>
   <td style="text-align:right;"> 60 </td>
   <td style="text-align:right;"> 89 </td>
   <td style="text-align:right;"> 106 </td>
   <td style="text-align:right;"> 72 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 44 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 85 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 95 </td>
   <td style="text-align:right;"> 110 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trnad-guc </td>
   <td style="text-align:right;"> 29 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 63 </td>
   <td style="text-align:right;"> 57 </td>
   <td style="text-align:right;"> 73 </td>
   <td style="text-align:right;"> 44 </td>
   <td style="text-align:right;"> 81 </td>
   <td style="text-align:right;"> 37 </td>
   <td style="text-align:right;"> 101 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 85 </td>
   <td style="text-align:right;"> 58 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 66 </td>
   <td style="text-align:right;"> 52 </td>
   <td style="text-align:right;"> 53 </td>
   <td style="text-align:right;"> 83 </td>
   <td style="text-align:right;"> 63 </td>
   <td style="text-align:right;"> 57 </td>
   <td style="text-align:right;"> 69 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 37 </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LOC111109073 </td>
   <td style="text-align:right;"> 52 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 41 </td>
   <td style="text-align:right;"> 44 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 52 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 42 </td>
   <td style="text-align:right;"> 122 </td>
   <td style="text-align:right;"> 63 </td>
   <td style="text-align:right;"> 79 </td>
   <td style="text-align:right;"> 74 </td>
   <td style="text-align:right;"> 45 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 92 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 68 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 25 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LOC111119351 </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 43 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 29 </td>
   <td style="text-align:right;"> 29 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 159 </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 42 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 46 </td>
   <td style="text-align:right;"> 44 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 50 </td>
   <td style="text-align:right;"> 31 </td>
   <td style="text-align:right;"> 75 </td>
  </tr>
</tbody>
</table>
  
**Thoughts** - Looks like most of these non-overlapping genes I have looked at so far are pseudo genes which don't have any known function and weren't in the transcriptome that was used for mapping by Salmon.  
  
### **Count Summaries**  

```r
# Column sums
st_col_sum <- colSums(st_c)
sa_col_sum <- colSums(sa_c)
# Column means
st_col_means <- colMeans(st_c)
sa_col_means <- colMeans(sa_c)

col_sums <- rbind(st_col_sum,sa_col_sum,st_col_means,sa_col_means)
row.names(col_sums) <- c("STAR - Count Sum","Salmon - Count Sum",
                         "STAR - Count Means","Salmon - Count Means")
colnames(col_sums) <- colnames(st_c)

kable(col_sums) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> RNA17005 </th>
   <th style="text-align:right;"> RNA17007 </th>
   <th style="text-align:right;"> RNA17013 </th>
   <th style="text-align:right;"> RNA17019 </th>
   <th style="text-align:right;"> RNA17069 </th>
   <th style="text-align:right;"> RNA17070 </th>
   <th style="text-align:right;"> RNA17072 </th>
   <th style="text-align:right;"> RNA17079 </th>
   <th style="text-align:right;"> RNA17090 </th>
   <th style="text-align:right;"> RNA17094 </th>
   <th style="text-align:right;"> RNA17099 </th>
   <th style="text-align:right;"> RNA17108 </th>
   <th style="text-align:right;"> RNA17122 </th>
   <th style="text-align:right;"> RNA17130 </th>
   <th style="text-align:right;"> RNA17142 </th>
   <th style="text-align:right;"> RNA17145 </th>
   <th style="text-align:right;"> RNA17162 </th>
   <th style="text-align:right;"> RNA17174 </th>
   <th style="text-align:right;"> RNA17176 </th>
   <th style="text-align:right;"> RNA17178 </th>
   <th style="text-align:right;"> RNA17181 </th>
   <th style="text-align:right;"> RNA17203 </th>
   <th style="text-align:right;"> RNA17211 </th>
   <th style="text-align:right;"> RNA17213 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> STAR - Count Sum </td>
   <td style="text-align:right;"> 185575.000000 </td>
   <td style="text-align:right;"> 212677.000000 </td>
   <td style="text-align:right;"> 215799.000000 </td>
   <td style="text-align:right;"> 207748.000000 </td>
   <td style="text-align:right;"> 250191.00000 </td>
   <td style="text-align:right;"> 202828.000000 </td>
   <td style="text-align:right;"> 186061.000000 </td>
   <td style="text-align:right;"> 197786.000000 </td>
   <td style="text-align:right;"> 205831.000000 </td>
   <td style="text-align:right;"> 184038.00000 </td>
   <td style="text-align:right;"> 206474.000000 </td>
   <td style="text-align:right;"> 183555.000000 </td>
   <td style="text-align:right;"> 216649.000000 </td>
   <td style="text-align:right;"> 189667.000000 </td>
   <td style="text-align:right;"> 170509.000000 </td>
   <td style="text-align:right;"> 173338.000000 </td>
   <td style="text-align:right;"> 184528.000000 </td>
   <td style="text-align:right;"> 262110.000000 </td>
   <td style="text-align:right;"> 196390.000000 </td>
   <td style="text-align:right;"> 156620.000000 </td>
   <td style="text-align:right;"> 156228.000000 </td>
   <td style="text-align:right;"> 192130.000000 </td>
   <td style="text-align:right;"> 186060.000000 </td>
   <td style="text-align:right;"> 219996.000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Salmon - Count Sum </td>
   <td style="text-align:right;"> 19574334.972000 </td>
   <td style="text-align:right;"> 20284503.975000 </td>
   <td style="text-align:right;"> 19283701.048000 </td>
   <td style="text-align:right;"> 21911886.995000 </td>
   <td style="text-align:right;"> 19555068.94300 </td>
   <td style="text-align:right;"> 19172671.971000 </td>
   <td style="text-align:right;"> 16892141.005000 </td>
   <td style="text-align:right;"> 20224865.017000 </td>
   <td style="text-align:right;"> 17484003.996000 </td>
   <td style="text-align:right;"> 16745087.05000 </td>
   <td style="text-align:right;"> 19283821.991000 </td>
   <td style="text-align:right;"> 16189340.990000 </td>
   <td style="text-align:right;"> 16612459.010000 </td>
   <td style="text-align:right;"> 17125343.932000 </td>
   <td style="text-align:right;"> 16145696.021000 </td>
   <td style="text-align:right;"> 18865380.972000 </td>
   <td style="text-align:right;"> 17772909.915000 </td>
   <td style="text-align:right;"> 20559587.991000 </td>
   <td style="text-align:right;"> 18737758.939000 </td>
   <td style="text-align:right;"> 16214359.967000 </td>
   <td style="text-align:right;"> 17797019.864000 </td>
   <td style="text-align:right;"> 15814310.993000 </td>
   <td style="text-align:right;"> 18239363.973000 </td>
   <td style="text-align:right;"> 16934144.937000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> STAR - Count Means </td>
   <td style="text-align:right;"> 4.757967 </td>
   <td style="text-align:right;"> 5.452837 </td>
   <td style="text-align:right;"> 5.532882 </td>
   <td style="text-align:right;"> 5.326462 </td>
   <td style="text-align:right;"> 6.41466 </td>
   <td style="text-align:right;"> 5.200318 </td>
   <td style="text-align:right;"> 4.770428 </td>
   <td style="text-align:right;"> 5.071046 </td>
   <td style="text-align:right;"> 5.277312 </td>
   <td style="text-align:right;"> 4.71856 </td>
   <td style="text-align:right;"> 5.293798 </td>
   <td style="text-align:right;"> 4.706176 </td>
   <td style="text-align:right;"> 5.554675 </td>
   <td style="text-align:right;"> 4.862882 </td>
   <td style="text-align:right;"> 4.371689 </td>
   <td style="text-align:right;"> 4.444222 </td>
   <td style="text-align:right;"> 4.731123 </td>
   <td style="text-align:right;"> 6.720252 </td>
   <td style="text-align:right;"> 5.035254 </td>
   <td style="text-align:right;"> 4.015588 </td>
   <td style="text-align:right;"> 4.005538 </td>
   <td style="text-align:right;"> 4.926031 </td>
   <td style="text-align:right;"> 4.770402 </td>
   <td style="text-align:right;"> 5.640489 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Salmon - Count Means </td>
   <td style="text-align:right;"> 512.860193 </td>
   <td style="text-align:right;"> 531.467078 </td>
   <td style="text-align:right;"> 505.245397 </td>
   <td style="text-align:right;"> 574.105562 </td>
   <td style="text-align:right;"> 512.35541 </td>
   <td style="text-align:right;"> 502.336363 </td>
   <td style="text-align:right;"> 442.584982 </td>
   <td style="text-align:right;"> 529.904499 </td>
   <td style="text-align:right;"> 458.092174 </td>
   <td style="text-align:right;"> 438.73207 </td>
   <td style="text-align:right;"> 505.248565 </td>
   <td style="text-align:right;"> 424.171169 </td>
   <td style="text-align:right;"> 435.257133 </td>
   <td style="text-align:right;"> 448.695049 </td>
   <td style="text-align:right;"> 423.027642 </td>
   <td style="text-align:right;"> 494.285141 </td>
   <td style="text-align:right;"> 465.661695 </td>
   <td style="text-align:right;"> 538.674457 </td>
   <td style="text-align:right;"> 490.941361 </td>
   <td style="text-align:right;"> 424.826682 </td>
   <td style="text-align:right;"> 466.293391 </td>
   <td style="text-align:right;"> 414.345141 </td>
   <td style="text-align:right;"> 477.883092 </td>
   <td style="text-align:right;"> 443.685512 </td>
  </tr>
</tbody>
</table>

### **Pearsons correlation between STAR and Salmon mappers (using only those sites that match**  

```r
st_c_match <- st_c[match(row.names(sa_c),row.names(st_c)),]
st_c_match$loc <- row.names(st_c_match)
sa_c_loc <- sa_c
sa_c_loc$loc <- row.names(sa_c)

cor_val <- rep(0,times=ncol(sa_c))
temp <- merge(st_c_match,sa_c_loc,by="loc")
for(i in 1:ncol(sa_c)){
  temp <- merge(st_c_match,sa_c_loc,by="loc")
  x <- temp[,c(i+1)]
  y <- temp[,c(i+25)]
  cor_val[i] <- cor(x,y)
}
```
  
Correlation by sample  

```r
names(cor_val) <- colnames(sa_c)
kable(cor_val) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> RNA17005 </td>
   <td style="text-align:right;"> 0.0875202 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17007 </td>
   <td style="text-align:right;"> 0.0768226 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17013 </td>
   <td style="text-align:right;"> 0.1270027 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17019 </td>
   <td style="text-align:right;"> 0.0757050 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17069 </td>
   <td style="text-align:right;"> 0.1203271 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17070 </td>
   <td style="text-align:right;"> 0.1101190 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17072 </td>
   <td style="text-align:right;"> 0.1197261 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17079 </td>
   <td style="text-align:right;"> 0.0861706 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17090 </td>
   <td style="text-align:right;"> 0.1081532 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17094 </td>
   <td style="text-align:right;"> 0.1212084 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17099 </td>
   <td style="text-align:right;"> 0.1032661 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17108 </td>
   <td style="text-align:right;"> 0.1206409 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17122 </td>
   <td style="text-align:right;"> 0.1190021 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17130 </td>
   <td style="text-align:right;"> 0.0930815 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17142 </td>
   <td style="text-align:right;"> 0.0795065 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17145 </td>
   <td style="text-align:right;"> 0.0877291 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17162 </td>
   <td style="text-align:right;"> 0.0667844 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17174 </td>
   <td style="text-align:right;"> 0.1386472 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17176 </td>
   <td style="text-align:right;"> 0.1060515 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17178 </td>
   <td style="text-align:right;"> 0.1110633 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17181 </td>
   <td style="text-align:right;"> 0.0924326 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17203 </td>
   <td style="text-align:right;"> 0.1338539 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17211 </td>
   <td style="text-align:right;"> 0.1221473 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNA17213 </td>
   <td style="text-align:right;"> 0.1241453 </td>
  </tr>
</tbody>
</table>
  
**Initial Thoughts** : The different mappers are not only leading to drastically different numbers of successfully mapped counts, but appear to lead to unique mapping to different transcripts that are only weakly correlated. In other words, Salmon leads not only to more reads being successfully mapped, but also leads to different loci with increased counts.  
  



