# Quick Analysis of Water Carbonate Chem Data in Multi-Dimensions (using PCA)

See for [Phenotype Analysis](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/markdown_files/Phenotype_Analysis/AE17_epfPhenotype.md) for complete analysis.

**Overview** : I used ```prcomp``` to decompose our three primary carbonate chem parameters ```EPF_pH```, ```EPF_DIC```, and ```EPF_Ca_Saturation```. We can see they are positively correlated in the first two PCs (not surprising, CA is determined using pH and DIC). Next I plotted our sample points along the first two PCs and color coded by either treatment or treatmentxtimepoint (shown below). We see that 2800xDay9 sticks out compared to the other treatmentxtime combinations. I am not sure we can used a PC to replace measured EPF_ph values given that we have so few points with all carbonated chemistry parameters measured.

![](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/notebook/img/PCA_waterChem_ExposureOnly_20190506.png)
