library(wasabi)
library(dplyr)

#Using Wasabi to convert Salmon files to ```.h5``` file format for sleuth

dir <- "/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files"
sfdirs <- filepaths("/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files/run20190610",list.files("/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files/run20190610"))
prepare_fish_for_sleuth(sfdirs)

model<-read.csv("/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files/metadata_cvirginica_rna_meta.txt", header=TRUE)
model$Treatment <- as.factor(model$treatment)
model$Time <- as.character(model$timepoint)
model$Time[model$Time == "3"] <- "09"
model$Time[model$Time == "6"] <- "80"
model$Time <- as.factor(model$Time)
model$Pop <- as.factor(model$population)
model$Lane <- as.factor(model$lane)
model$SFV <-  interaction(model$Time,model$Treatment) # Creates single factor variable for combination of time and treatment
samp_names <- substr(model$sample_name,4,10)
model <- dplyr::mutate(model, path = file.path(dir,'/run20190610_h5/',samp_names,'_abundance.h5',fsep = ""))

metadata <- dplyr::select(model,c("sample_name","Time","Treatment","SFV","Pop","Lane","path"))
colnames(metadata)[1] <- "sample"

trans <- readRDS("/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files/NCBI_transcriptome/transcriptome_table.RData")

ttg <- dplyr::select(trans,c("fullID","location","product"))
colnames(ttg) <- c("target_id","ens_gene","product")

so <- sleuth_prep(metadata, target_mapping = ttg,aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE)
saveRDS(so,"/shared_lab/20180226_RNAseq_2017OAExp/RNA/salmon_files/20190610_sluethobject_fullmodelparameters.RData")
