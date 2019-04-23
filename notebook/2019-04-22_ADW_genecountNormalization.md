# Revised gene count normalization script

Update [03B_CV17_RNA_countNormalization]()

Script updated to take ```scenario1``` rdata object from the most recent verions of the [filtering script]() (contains genes with less than 50 counts (```CT < 50```) or genes with the majority of counts in a single gene (```PMAX < 0.5```)), then normalizes and data in a two step process.
* Step one: uses ```TMM - Normalization``` implemented using ```calcNormFactors``` in the ```edgeR``` package to normalize data by accoutning for between sample variability.
* Step two: log normalizes data with ```voomWithQualityWeights()``` function in the ```limma``` package. The choice was also made to include tank as a blocking variable in this function.

The script creates at final normalized .RData file (e.g. ```scenarioX_NormalizedVoom_Rdata```) in the ```input_file/RNA``` folder directory.
* The current version only takes the data from ```scenario1```

