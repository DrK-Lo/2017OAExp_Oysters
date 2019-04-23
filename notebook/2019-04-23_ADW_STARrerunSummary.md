# Reran STAR mapping program

To double check the performance of the STAR mapping program and the original read alignment step I reran STAR using the same settings Brett used when he originally ran the program. I had to rewrite a few of the scripts to accomplish this, but the source data (trimmed reads and genome files) remained the same between runs. Furthermore, the arguements remained unchanged. I performed the same two step process storing the data in a new folder called ```/STAR_output_v2``` in the ```RNA``` directory```. Within the folder there are the two subfolder (```m2``` and ```m3```) for the two different passes.

The second pass is still finishing, but the initial observation is that resulting mapping reads have output files that are potentiall close but NOT identical. In particular, two observations that appear systemic across all individuals mapped thus far and between both passes:
* 1 - The most recent run starts with about ~10 million more reads 
* 2 - There are more overall reads mapped to single loci, but the percentage itself is lower (consequently the multi-map percentage is a little higher)

**Example - Sample RNA17005**

Side by side comparison of the two runs. (Most recent run on the left)

M2
![](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/figures/m2_STAR_outputLog_comparison.png)

M3
![](https://github.com/epigeneticstoocean/2017OAExp_Oysters/blob/master/figures/m3_STAR_outputLog_comparison.png)

**Preliminary Conclusions**
Unclear why we are seeing this difference given that the data and arguements are the same between runs. I will be following up with brett to see if there is any reasong he can think of for these differences.
