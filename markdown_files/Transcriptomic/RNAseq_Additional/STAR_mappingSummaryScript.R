# Script for taking STAR mapping summary file and transforming it into a readable table

# Code in bash used to generate STAR mapping summary file.
# grep -r 'Number of input reads\|Uniquely mapped reads number\|Uniquely mapped reads %\|Number of reads mapped to multiple loci\|% of reads mapped to multiple loci' *.final.out > STAR_simpleMappingSummaryAllSamples.txt


star_df <- read.delim("/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/extra/STAR_simpleMappingSummaryAllSamples.txt",header=FALSE,stringsAsFactors = FALSE)
star_list <- strsplit(star_file[,1],":")
star_mat <- do.call(rbind,star_df)
star_labels <- unique(substr(star_mat[,1],4,8))
star_df_mat <- matrix(star_df[,2],ncol=5,byrow = TRUE)
colnames(star_df_mat) <- c("Input_Reads","Unique_Reads","Unique_Reads_percent","Multi_Reads","Multi_Reads_percent")
star_final_df <- data.frame(ID=star_labels,star_df_mat)

write.csv(star_final_df,"/home/downeyam/Github/2017OAExp_Oysters/input_files/RNA/extra/STAR_mappingSummaryFinalTable.csv")
