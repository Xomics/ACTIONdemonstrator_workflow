#!/usr/bin/Rscript

################################################################################
## Title:         heatmap_missingness.R
## Description:   Create heatmaps that visualize NA values in all datasets
## Author:        Casper de Visser
## Date created:  2022-07-02
## Email:         casper.devisser@radboudumc.nl
################################################################################
## Notes:
## 
################################################################################



args = commandArgs(trailingOnly=TRUE)

number_of_arguments <- length(args)
number_of_omics <- number_of_arguments-1
output_dir = args[number_of_arguments]


# Read all csv files
omicsData <- lapply(args[1:number_of_omics], function(x) {read.csv(x, row.names=1)})



# Functions

#'  Subset specific columns of dataframe or matrix
#'  
#'  @param data input data
#'  @param cols cooridates of columns that are subsetted
#'  @return subsetted data
subset_col_string <- function(data, string) {
  subset <- data[, grepl(string, colnames(data))]
  return(subset)
}

#'Function to remove ACTIONBB3/ACTIONBB23 from headers
#'
#' @param df Metabolomics df of one assay/platform
#' @param assay The name of the assay/platform as as string
#' @return df of one assay with colummn names without ACTIONBB3 + assay name
remove_ACTION_in_headers <- function(df, study, assay){
  new_colnames <- lapply(colnames(df), function(x) gsub(paste0(study, '_', assay, '.'), '', x))
  colnames(df) <- new_colnames
  return(df)
}


#' Pre-processing of heatmaps of missingness of metabolomics data
#'
#' @param data metabolomics_values
#' @return list of metabolomics subset(amines, OA, steroids, Bio and Dipstick) with headers rewritten as solely the metabolites
heatmap_missingness_preprocessing_metabolomics <- function(data) {
  #Get list of all metabolomics platforms
  platform_subsets <- lapply(METABOLOMICS_PLATFORMS, function(x) subset_col_string(data, x))
  
  #Remove 'ACTIONBB3' from headers, leaving only the metabolites
  platform_subsets_new_headers <- lapply(seq(1, 3), function(x) remove_ACTION_in_headers(platform_subsets[[x]], 'ACTIONBB3', METABOLOMICS_PLATFORMS[[x]]))
  
  #Subset only ACTIONBB23 columns of Bio and Dipstick
  Bio_ACTIONBB23 <- subset_col_string(platform_subsets[[4]], 'ACTIONBB23')
  DS_ACTIONBB23 <- subset_col_string(platform_subsets[[5]], 'ACTIONBB23')
  
  #Remove 'ACTIONBB23' from headers, leaving only the metabolites
  Bio_ACTIONBB23_new_headers <- remove_ACTION_in_headers(Bio_ACTIONBB23, 'ACTIONBB23', METABOLOMICS_PLATFORMS[[4]])
  DS_ACTIONBB23_new_heeaders <- remove_ACTION_in_headers(DS_ACTIONBB23, 'ACTIONBB23', METABOLOMICS_PLATFORMS[[5]])
  
  list_of_subsets <- list(platform_subsets_new_headers[[1]], 
                          platform_subsets_new_headers[[2]], 
                          platform_subsets_new_headers[[3]],
                          Bio_ACTIONBB23_new_headers,
                          DS_ACTIONBB23_new_heeaders)
  return(list_of_subsets)
}



#' Heatmap of missingness data
#' 
#' @param Data Dataframe of which missingness will be represented in heatmap
#' @return Heatmap of the missingness
heatmap_missing <- function(data, title) {
  library(ggplot2) 
  heatmap <- visdat::vis_miss(data, show_perc_col = FALSE, warn_large_data = FALSE) + 
                ggtitle ('Missing values', title) +
                theme(plot.title = element_text(size = 16, hjust = 0.5)) +
                if (ncol(data) > 50) {
                  theme(axis.text.x = element_blank()
                   )}
  print(heatmap) 
  return(heatmap)
}


path_splits <- lapply(args[1:number_of_omics], function(x) {strsplit(x, "/")})
file_names <- lapply(path_splits, function(x) tail(x[[1]], n=1))

p <-  lapply(seq(1:number_of_omics), function(x) heatmap_missing(omicsData[[x]], file_names[[x]]))



pdf(output_dir)
p
dev.off()
