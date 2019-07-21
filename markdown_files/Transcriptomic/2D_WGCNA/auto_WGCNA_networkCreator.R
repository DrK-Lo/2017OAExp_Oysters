
setwd("/shared_lab/20180226_RNAseq_2017OAExp/RNA/WGCNA");
library(WGCNA)
library(tximport)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. 
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads(nThreads = 40)
# Load the data saved in the first part
lnames = load(file = "20180513_salmonrun_consensus.RData");
# The variable lnames contains the names of loaded variables.
#lnames

writeLines("Starting net function....")
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 15,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "20180513_SalmonRunTOM", 
                       verbose = 3)
writeLines("Net function complete, storing results....")

moduleLabels = net$colors;
moduleColors = labels2colors(net$colors);
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "20180513_SalmonRun_network_auto.RData");