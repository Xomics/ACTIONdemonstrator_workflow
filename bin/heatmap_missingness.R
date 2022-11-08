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
