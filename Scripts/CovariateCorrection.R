###################################################################################
#
# 05/08/2022
#
# Project: x-omics ACTION Demonstrator - FAIR multi-omics workflow
# script: Example code correct omics data for covariates
#
###################################################################################

rm(list=ls())
#r-packages
library(nlme) #R-package for (non)linear mixed-effects models
library(foreach) #R-package I prefer over straight forward for loops 
library(doMC) #allows for parallelization of foreach functions
registerDoMC(X) #set number of cpu's to use in foreach functions by replacing X with desired/available number

###################################################################################
#
# Example code correction 
# without accounting for clustering in the data
# i.e., no correction for relatedness/zygosity
#
# assumptions:
# 1) omics data are contained in dataframe called omics
# 2) covariates are contained in dataframe called covs
# 3) N rows omics & covs dataframes are equal + order of rows are also equal
# 4) names of omics data to correct are contained in vector called omics.traits
#
###################################################################################

args = commandArgs(trailingOnly=TRUE)

omics_path = args[1]
covariates_path = args[2]
omics_output = args[3]

omics <- read.csv(omics_path, row.names=1)
covs <- read.csv(covariates_path, row.names=1)



# Remove 'X' from omics columns names #TODO: this is only neaded for real data, thinkg of long term solution
if (colnames(omics)[1] != "XOE1") {
   colnames(omics) <- gsub("X", "", colnames(omics)) #TODO: Hardcoded, put in other epigenomics pre-processing script?
}

# Transforms rows and columns
omics <- t(omics)

#for each omics trait calculate the residual after correction for covariates (results combine into matrix)
omics_corrected <- resid(lm(omics ~ covs$Sample_Plate + covs$Array_rownum + covs$Epi + covs$Fib + covs$B + covs$NK + covs$CD4T + covs$CD8T + covs$Mono + covs$Neutro+ covs$Eosino , na.action=na.exclude)) # obtain residuals from linear model incl. covariates (used sex + age as example covariates here)

#matrix to dataframe
omics_corrected <- as.data.frame(omics_corrected) #matrix2df

# Transform rows and columns back
omics_corrected <- t(omics_corrected)

#write to csv
write.csv(omics_corrected, omics_output)
