#!/usr/bin/Rscript

################################################################################
## Title:         concatenate_MAF.R
## Description:   Concatenate metabolomics platforms into single file
## Author:        Casper de Visser
## Date created:  2022-05-23
## Email:         casper.devisser@radboudumc.nl
################################################################################
## Notes:
## 
################################################################################


args = commandArgs(trailingOnly=TRUE)

mtblmcs_amines_path = args[1]
mtblmcs_oa_path = args[2]
mtblmcs_steroids_path = args[3]

output_dir = args[4]

amines <- read.csv(mtblmcs_amines_path, header = TRUE)
oa <- read.csv(mtblmcs_oa_path, header = TRUE)
steroids <- read.csv(mtblmcs_steroids_path, header = TRUE)


list1 <- c(nrow(amines), nrow(oa), nrow(steroids))

# Check if number of rows (samples) is equal. Find common elements if not
if (max(list1) != min(list1)) {
  common_IDs <- intersect(intersect(amines[,1], oa[,1]), steroids[,1])
  
  amines <- subset(amines, subset = amines[,1] %in% common_IDs)
  oa <- subset(oa, subset = oa[,1] %in% common_IDs)
  steroids <- subset(steroids, subset = steroids[,1] %in% common_IDs)
  
  binded_df <- cbind(amines, oa[,2:ncol(oa)], steroids[,2:ncol(steroids)])
} else {
  
  binded_df <- cbind(amines, oa[,2:ncol(oa)], steroids[,2:ncol(steroids)])
  
}


write.csv(binded_df, file = output_dir, row.names = FALSE )
