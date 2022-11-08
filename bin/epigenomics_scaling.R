#!/usr/bin/Rscript

################################################################################
## Title:         epigenomics_scaling.R
## Description:   Normalization of epigenomics data
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


#'  Normalization data 
#'  
#'  @param data The input data that requires normalization
#'  @return Data matrix normalized by mean centering 
normalization <- function(data) { 
  data_mean_centered <- scale(data, scale = TRUE, center = TRUE)
  return(data_mean_centered)
}


# Transform rows and columns
df <- t(df)

# Normalization of data for PCA
epigenetics_normalized <- normalization(df)

# Transform rows and columns back
epigenetics_normalized <- t(epigenetics_normalized)
  

# Save df to .csv
write.csv(epigenetics_normalized, file = output_dir, row.names = TRUE)
