# Analysis of EPF Fluid Data for 2017 OA Exposure Experiment

**Overview** : Following the plans outlined in the [previous entry](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/notebook/20190428_EPFanalysisPlan.md), I reran the EPF data analysis testing for a significant relationship between EPF fluid chemistry (mostly pH, some other carbonate chem parameters when available) and time or treatment. I performed this on three different subsetted versions of the data ( a complete workup of the data can be found [here](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/Phenotype_Analysis/AE17_epfPhenotype.md)). Below I will focus on two subsets:
* First :  Using all data (6 timepoints, 3 treatment levels, 104 observations)
* Second :  Only data for sequenced individuals (2 timepoints, 2 treatments, 24 observations)

### **Complete Dataset analysis**

* Response Variable : EPF pH
* Explanatory Variables : Time (24hr,36hr,9days,22days,50days,80days) and Treatment (400,900,2800)
* Random Effects : Population (3 'pops'), Shelfs (6 shelves), Tanks (18 tanks; nested in shelves)

![](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/notebook/img/fullEPFpH_boxplot_20190502.png)

**Full Model(s)** (show in lmer format)

* Two-way ANOVA
  * EPF pH ~ Time(factor) * Treatment + (1|Population) + (1|Shelf/Test)
* ANCOVA
  * EPF pH ~ Time(cont) * Treatment + (1|Population) + (1|Shelf/Test)
  
**Final Models**

* Two-way ANOVA
  * EPF pH ~ Time(factor) * Treatment + (1|Test:Shelf)
* ANCOVA
  * EPF pH ~ Treatment + (1|Test:Shelf)
  
**Final Analysis**  

Summary table for the best Two-way ANOVA table. This is the only one presented given that the ANCOVA actually ideally simplified into a model nested in the model selected for the Two-way ANOVA. 

![](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/notebook/img/full_EPFpH_modelTable_output_20190502.png)
  
#### **Brief Description of Process**

*IN PROGRESS, but see full .md for complete details
  
### **Sequenced Sample Dataset analysis**

* Response Variable : EPF pH
* Explanatory Variables : Time (9days, 80days) and Treatment (400,2800)
* Random Effects : Population (3 'pops'), Shelfs (6 shelves), Tanks (18 tanks; nested in shelves)

**Bar plot of EPF pH with error bars and final significance levels**
![](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/notebook/img/sequencec_EPFpH_barplot_20190502.png)

**Full Model(s)** (show in lmer format)

* Two-way ANOVA
  * EPF pH ~ Time(factor) * Treatment + (1|Population) + (1|Shelf/Test)
  
**Final Models**

* Two-way ANOVA
  * EPF pH ~ Time(factor) * Treatment
  
**Final Analysis**  

Summary table for the best Two-way ANOVA table.

![](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/notebook/img/sequence_EPFpH_table_modelOutput_20190502.png)

**NOTE**: This analysis is for a dataset with three outliers removed (see full markdown for details)
  
#### **Brief Description of Process**

*IN PROGRESS, but see full .md for complete details
  
