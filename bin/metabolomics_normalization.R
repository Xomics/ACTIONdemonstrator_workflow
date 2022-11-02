#!/usr/bin/Rscript

################################################################################
## Title:         metabolomics_normalization.R
## Description:   Normalize metabolomics data
## Author:        Casper de Visser, René Pool
## Date created:  2022-02-02
## Email:         casper.devisser@radboudumc.nl
################################################################################
## Notes:
## 
################################################################################


args = commandArgs(trailingOnly=TRUE)

mtblmcs_filtered_path = args[1]
mtblcs_values_path = args[2]

output_dir = args[3]

normalization_method = args[4]


# Read metabolomics filtered data
metabolomics_NA_filtered <- read.csv(mtblmcs_filtered_path)


# Read metabolomcics data (for bio and DS values)
load_file <- function(infile) {
	require(tools)
	if (file_ext(infile) == "csv") {
		df <- read.csv(infile)
	} else if (file_ext(infile) == "tsv") {
		df <- read.csv(infile, sep = "\t")
	}
	return(df)
}

mtblcs_values <- load_file(mtblcs_values_path)



#' Performs metabolomics normalization based on measured creatinine levels
#' 
#' @param metabolomics_values_df dataframe containing the metabolomic data
#' @return CRNormalizedMtblmcsDF dataframe where the values have been normalized
cr_normalize_metabolomics_values <- function(metabolomics_values_df){
  
  # Normalizes the intensities of the metabolites by the creatinine levels measured at GBS
  #
  # The goal of this operation is correcting for the hydration state of the subjects, enabling
  # the study of the relation between metabolic variation between individuals and a phenotype 
  # (in our case aggression)
  #
  # We will use variable 1 / ACTIONBB23_bio_creatinine as CR correction factor
  # 
  # FCR = 1/ACTIONBB23_bio_creatinine
  #
  # I_corrected = I_reported * FCR.
  # 
  # See also "../../Documents/Miller 2004 SG and creatinine urine correction.pdf" (https://doi.org/10.1373/clinchem.2004.032292).
  
  # Remove metabolites that were filtered out in pre-processing
  mtblcs_values <- subset(mtblcs_values, mtblcs_values[,1] %in% metabolomics_NA_filtered$XOmicsmetaboID)
  creatinine_data <- mtblcs_values$ACTIONBB23_bio_creatinine
  
  
  MtblmcsDF <- metabolomics_values_df
  Columns <- names(MtblmcsDF[,2:ncol(MtblmcsDF)])
  #MtblmcsDF$XOmicsmetaboID <- row.names(MtblmcsDF)
  MtblmcsDF <- MtblmcsDF[,
                         c("XOmicsmetaboID",
                           Columns)]
  remove(Columns)
  MtblmcsDF$FCR <- 1.0 / creatinine_data
  
  VariablesToBeNormalized <- names(MtblmcsDF[,2:ncol(MtblmcsDF)])
  
  # Code is redundant, when using the MAF files as input (which already separate the plaforms)
  # Filter <- grepl(pattern = "_amines",
  #                 x = VariablesToBeNormalized)
  # Filter <- Filter | grepl(pattern = "_OA",
  #                          x = VariablesToBeNormalized)
  # Filter <- Filter | grepl(pattern = "_steroids",
  #                          x = VariablesToBeNormalized)
  #VariablesToBeNormalized <- VariablesToBeNormalized[Filter]
  #remove(Filter)
  
  # CRNormalizedMtblmcsDF <- MtblmcsDF[,
  #                                    c("XOmicsmetaboID",
  #                                      VariablesToBeNormalized,
  #                                      "FCR")]
  
  CRNormalizedMtblmcsDF <- MtblmcsDF
  
  for(V in VariablesToBeNormalized){
    CRNormalizedMtblmcsDF[,
                          V] <- 
      CRNormalizedMtblmcsDF[,
                            V] * CRNormalizedMtblmcsDF$FCR
  }
  CRNormalizedMtblmcsDF$FCR <- NULL

  return(CRNormalizedMtblmcsDF)
}



#' Performs metabolomics normalization based on specific gravity
#' 
#' @param metabolomics_values_df dataframe containing the metabolomic data
#' @return SGNormalizedMtblmcsDF dataframe where the values have been normalized
sg_normalize_metabolomics_values <- function(metabolomics_values_df){
  
  # Normalizes the intensities of the metabolites by the specific gravity levels measured at GBS
  #
  # The goal of this operation is correcting for the hydration state of the subjects, enabling
  # the study of the relation between metabolic variation between individuals and a phenotype 
  # (in our case aggression)
  #
  # We will use variable ACTIONBB23_bio_Density for creating an SG correction factor
  # 
  # FSG = (SG_target - 1) / (SG_sample - 1),
  # where SG_target = median(ACTIONBB23_bio_Density) and SG_sample = ACTIONBB23_bio_Density.
  #
  # I_corrected = I_reported * FSG.
  # 
  # See also "../../Documents/Miller 2004 SG and creatinine urine correction.pdf" (https://doi.org/10.1373/clinchem.2004.032292).
  
  # Remove metabolites that were filtered out in pre-processing
  mtblcs_values <- subset(mtblcs_values, mtblcs_values[,1] %in% metabolomics_NA_filtered$XOmicsmetaboID)
  density_data <- mtblcs_values$ACTIONBB23_bio_Density
  
  MtblmcsDF <- metabolomics_values_df
  Columns <- names(MtblmcsDF[,2:ncol(MtblmcsDF)])
  MtblmcsDF$XOmicsmetaboID <- row.names(MtblmcsDF)
  MtblmcsDF <- MtblmcsDF[,
                         c("XOmicsmetaboID",
                           Columns)]
  remove(Columns)
  

  Filter= (density_data == 1.0)
  # table(Filter)
  # There is 1 sample with density == 1
  #   => set this value to NA
  density_data[Filter] <- NA
  remove(Filter)

  
  MtblmcsDF$FSG <- median(x = density_data,
                          na.rm = TRUE) - 1.0
  MtblmcsDF$FSG <- MtblmcsDF$FSG / (density_data - 1.0)
  
  VariablesToBeNormalized <- names(MtblmcsDF[,2:ncol(MtblmcsDF)])
  
  # Code is redundant, when using the MAF files as input (which already separate the plaforms)
  # Filter <- grepl(pattern = "_amines",
  #                 x = VariablesToBeNormalized)
  # Filter <- Filter | grepl(pattern = "_OA",
  #                          x = VariablesToBeNormalized)
  # Filter <- Filter | grepl(pattern = "_steroids",
  #                          x = VariablesToBeNormalized)
  # VariablesToBeNormalized <- VariablesToBeNormalized[Filter]
  # remove(Filter)
  
  # SGNormalizedMtblmcsDF <- MtblmcsDF[,
  #                                    c("XOmicsmetaboID",
  #                                      VariablesToBeNormalized,
  #                                      "FSG")]
  
  SGNormalizedMtblmcsDF <- MtblmcsDF
  
  for(V in VariablesToBeNormalized){
    SGNormalizedMtblmcsDF[,
                          V] <- 
      SGNormalizedMtblmcsDF[,
                            V] * SGNormalizedMtblmcsDF$FSG
  }
  SGNormalizedMtblmcsDF$FSG <- NULL

  return(SGNormalizedMtblmcsDF)
}





# Normalize data, method selected with Nextflow param
if (normalization_method == 'cr'){
  metabolomics_values_normalized <- cr_normalize_metabolomics_values(metabolomics_NA_filtered)

 }

if (normalization_method == 'sg'){
  metabolomics_values_normalized <- sg_normalize_metabolomics_values(metabolomics_NA_filtered)

} 


# Write normalized files to .csv

write.csv(metabolomics_values_normalized, file = output_dir, row.names = FALSE)
