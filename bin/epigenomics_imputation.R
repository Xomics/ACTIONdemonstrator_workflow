#!/usr/bin/Rscript

################################################################################
## Title:         epigenomics_imputation.R
## Description:   Imputation and normalization of epigenomics data
## Author:        Casper de Visser, Jenny van Dongen
## Date created:  2022-03-02
## Email:         casper.devisser@radboudumc.nl
################################################################################
## Notes:
## 
################################################################################


#Read in data

args = commandArgs(trailingOnly=TRUE)

epigenomics_values_path = args[1]
output_dir = args[2]

df <- read.csv(epigenomics_values_path, row.names = 1)


#'  Imputation epigenetics data
#' 
#' @param unprocessed_data Unprocessed dataframe for PCA
#' @return Data matrix with only numerical and NA's imputated with median
data_imputation_median <- function(unprocessed_data) {
  nums <- unlist(lapply(unprocessed_data, is.numeric))
  data_num <- as.matrix(unprocessed_data[,nums])
  data_num[is.na(data_num)] <- median(data_num, na.rm = TRUE)
  return(data_num)
}



# Data imputation epigenetics
epigenetics_imputed <- data_imputation_median(df)

 
  

# Save df to .csv
write.csv(epigenetics_imputed, file = output_dir, row.names = TRUE)
