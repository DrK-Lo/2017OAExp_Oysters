# DEG Analysis using RSEM estimated genes and Limma voom method


Output from RSEM analysis was aggregated into a single gene count matrix (columns = samples, rows = gene LOC).

**Data**
* RSEM count matrix (R object)
  * `~/2017OAExp_Oysters/input_files/RNA/STAR_mapping/RSEM_output/RSEM_gene_Summary.Rdata`
* Meta datasheet for each sample
  * `/2017OAExp_Oysters/input_files/meta/metadata_20190718.RData`
* Custom gene annotation file created for use with tximport and based on gtf file (with only Gnomon based features) used for the original STAR index and mapping
  * `~/2017OAExp_Oysters/input_files/RNA/references/STAR_gnomon_tximportGeneFile.RData`

**Libraries**
* `limma`
* `edgeR`

### Reading in data
```
#### RSEM counts ####
RSEM <-  readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/STAR_mapping/RSEM_output/RSEM_gene_Summary.Rdata")
# Separate out RSEM counts and rename rows with LOC ID
rsem_c <- RSEM$Count

#### Transcript File ####
tran <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/references/STAR_gnomon_tximportGeneFile.RData")
# Gene File
gene <- tran[!duplicated(tran$GENEID),]

#### Meta Data ####
meta <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/meta/metadata_20190718.RData")
meta$sampl_nameSimple <- substr(meta$sample_name,start = 4,stop=9)

#### Data Manipulation ####
# Separate out RSEM counts and rename rows with LOC ID
loc <- gene$GENEID[match(row.names(rsem_c),gene$gene_id)]
row.names(rsem_c) <- loc
rsem_order <- rsem_c[order(rownames(rsem_c)),] # order by LOC id

#### Create DGEList obj (using EdgeR package)
rsem <- DGEList(rsem_order) 
```

### Step One : Filtering and Normalization*
* Filtering removed genes with limited to no expression
  * CRITERIA: removed individuals that have less than 1 count per million (CPM) in less than have the total number of samples.
* Normalization was done with the `limma` function : `calcNormFactors()`

**Code**
```
keep_rsem <- rowSums(cpm(rsem)>1) >= (0.5 * 24)
rsem_red <- rsem[keep_rsem, ]
rsem_red <- calcNormFactors(rsem_red)
```

### DEG Analysis using basic planned comparisons

**NOTE**: Details for performing planned comparisons in `limma` can be found pg ~48 of the [user guide](https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)

**Primary Questions**
1) which genes are responsive to OA after short exposure
2) which genes are responsive to OA after long exposure
3) which genes are consistently responsive to OA (regardless of exposure duration)
4) which genes are sensitive to the duration in our experimental system under control conditions

**Planned comparisons**
1) 400_Day09 - 2800_Day09
2) 400_Day80 - 2800_Day80
3) 400 - 2800
4  400_Day09 - 400_Day80

Standard Design Model in `limma`
```
design <- model.matrix(~Treatment*Time)
```

Expanded Design
| coefficient | Comparison | Interpretation |
|:-----------:|:----------:|:--------------:|
| Intercept | 400.09 | Baseline level of unexposed at Day 09 |
| Treatment09 | 400.09-2800.09 | Difference between exposures at Day 09 |
| Time400 | 400.09-400.80 | Effect of duration (9 vs 80) in unexposed |
| Treatment09:Time400 | (2800.80 - 2800.09)-(400.80-400.09) | Interaction |

**Run the model (count matrix, model design)**
```
rsem_fit <- lmFit(rsem_out, design)
```

*NOTE*: Not all planned comparisons are achieved when running the standard model.

Using contrasts to examine all specific comparisons.


