#!/usr/bin/Rscript

################################################################################
## Title:         metabolomics_scaling.R
## Description:   Pareto scaling of metabolomics data
## Author:        Casper de Visser
## Date created:  2022-10-19
## Email:         casper.devisser@radboudumc.nl
################################################################################
## Notes:
## 
################################################################################


#Read in data

args = commandArgs(trailingOnly=TRUE)

metabolomics_values_path = args[1]
output_dir = args[2]


# Read data
df <- read.csv(metabolomics_values_path)



#'  Pareto scaling of data 
#'  Function modified from the RFmarkerDetector package
#'  
#'  @param data The input data that needs to be scaled
#'  @return Data matrix scaled
paretoscale <- function(data) {
   x <- data
   x.centered <- apply(x, 2, function(x) x - mean(x, na.rm = T))
   x.sc <- apply(x.centered, 2, function(x) x/sqrt(sd(x, na.rm = T)))
   data.frame(x.sc)
}


# Subset numerical data
df_num <- df[,2:ncol(df)]

# Scale data
metabolomics_scaled <- paretoscale(df_num)

# Add samples names column
merged_dfs <- cbind(df[,1], metabolomics_scaled)
colnames(merged_dfs) <- colnames(df)

# Save df to .csv
write.csv(merged_dfs, file = output_dir, row.names = FALSE)
