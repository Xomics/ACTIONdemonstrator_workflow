#!/usr/bin/Rscript

################################################################################
## Title:         sort_cols_sd.R
## Description:   Sort columns on standard deviation, then take subset columns (10%) with highest variance
## Author:        Casper de Visser
## Date created:  2022-28-03
## Email:         casper.devisser@radboudumc.nl
################################################################################
## Notes:
## 
################################################################################


#Read in data

args = commandArgs(trailingOnly=TRUE)

epigenomics_values_path = args[1]
feature_cutoff = args[2]
output_dir = args[3]


df <- read.csv(epigenomics_values_path, row.names = 1)


#' Get standard deviation of column with NA
#'
#' @param column of numericals
#' @return standard deviation
get_sd <- function(column) {
	SD <- sd(column, na.rm = TRUE)
	return(SD)
}


#' Subset dataframe based on SD of columns
#' 
#' @param dataframe any dataframe (with numericals)
#' @return Same dataframe, only containing columns that belong to top 10% with regard to standard deviation
subset_cols_sd <- function(df) {
  standard_deviations <- apply(df, 2, get_sd)
  standard_deviations_sorted <- sort(standard_deviations, decreasing = TRUE)
  ten_percent  <- ncol(df) / as.numeric(feature_cutoff)
  names <- names(standard_deviations_sorted[1:ten_percent])
  
  data_subset <- df[,names]

  return(data_subset)
}


# Transform row and columns
df <- t(df)

# Subset features CpG sites with highest variance
data_subset <- subset_cols_sd(df)

# Transform rows and columns back
data_subset <- t(data_subset)

# Save df to .csv
write.csv(data_subset, file = output_dir, row.names = TRUE)
