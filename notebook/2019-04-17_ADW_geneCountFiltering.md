# Updated filtering script

Updated [03A_CV17_RNA_countFilteringandAnalysis.Rmd](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/STAR_pipeline/03A_CV17_RNA_countFilteringandAnalysis.md)

Script looks at different filtering approaches to remove uninformative or problematic genes.

It includes:
 * A number of diagnostic plots (looking specifically at lane bias, proportion of counts explained by a single individual, number of samples with counts per gene, and a few other measures of gene diversity.)
 * Filtered the full geneCount and stored new matrix for downstream analysis.
  * Currently includes a single ```scenario1``` which includes filtering by a minimum count ```CT = 50```` and a maximum proportion of counts ins a single individuals ```PMAX = 0.5```
 *```Scenario1``` stored as an Rdata object in the ```inpute_files/RNA``` folder. 

