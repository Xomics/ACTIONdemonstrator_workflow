#!/usr/bin/Rscript

################################################################################
## Title:         pca.R
## Description:   Normalize metabolomics data
## Author:        Casper de Visser, Jenny van Dongen
## Date created:  2022-02-02
## Email:         casper.devisser@radboudumc.nl
################################################################################
## Notes:
## 
################################################################################

args = commandArgs(trailingOnly=TRUE)

mtblmcs_values_path = args[1]
epigenomics_values_path = args[2]
epigenomics_meta_path = args[3]
phenotypes_set2_path = args[4]
ids_path = args[5]

output_dir = args[6]


# Load data

# Read metabolomics preprocessed data
metabolomics_values_normalized <- read.csv(mtblmcs_values_path, header = TRUE, row.names = 1)


# Read epigenomics preprocessed data
epigenetics_preprocessed <- read.csv(epigenomics_values_path, row.names = 1, header = TRUE)

# Epigenomics metadata
epigenetics_meta <- read.csv(epigenomics_meta_path, row.names = 1, header = TRUE)


# Read phenotype set2
phenotypes_set2 <- read.csv(phenotypes_set2_path)
phenotypes_set2 <- data.frame(phenotypes_set2, row.names=1)

# Read IDs fle
IDs <- read.table(ids_path, header = TRUE, sep = ",")


# Libraries
library(ggplot2)


# Functions

#'  Filter samples (rows) with missing values
#' 
#' @param df Dataframe with missing values
#' @return Data matrix with only numericals and no missing values
filter_samples_missing_values <- function(df) {
  row_names_with_NA <- rownames(df[rowSums(is.na(df)) > 0,])
  data_without_NA <- df[!row.names(df) %in% row_names_with_NA,]
  return(data_without_NA)  
}



#' Imputation of missing values (using median)
#'
#' @param df Dataframe with missing values (samples as columns, features as rows)
#' @return df without missing values
data_imputation_median <- function(df) {
  nums <- unlist(lapply(df, is.numeric))
  data_num <- as.matrix(df[, nums])
  data_num[is.na(data_num)] <- median(data_num, na.rm = TRUE)
  return(data_num)
}



#' Perform PCA
#' 
#' @param preprocessed_data Data matrix that is preprocessed and normalized
#' @param output_dir The directory where pca .RData is stored
#' @param name_file Name of the .RData file
#' @return The PCA
perform_PCA <- function(preprocessed_data) {
  pc <- prcomp(preprocessed_data, scale. = TRUE)
  return(pc)
}


#'  Create a screeplot
#'  
#'  @param pc PCA object as input
#'  @return The screepot object 
make_screeplot <- function(pc) {
  
  # create summary of pca
  vars <- pc$sdev^2
  vars <- vars/sum(vars)
  df <- as.data.frame(cbind('SD' = pc$sdev, 'Prop' = vars, 'cumulative'  = cumsum(vars)))
  df$PC <- seq(nrow(df))
  df$Prop <- df$Prop * 100
  df$cumulative <- df$cumulative * 100
  df <- df[c('PC', 'SD', 'Prop', 'cumulative')]
  df <- as.data.frame(lapply(df, unlist))
  
  # Get number of PCs that explain 80%
  eighty_percent <-  min(which(df$cumulative > 80))
  scale <- round(eighty_percent/25)
  if (scale == 0) {
    scale <- 1
  }
  
  # Barplot
  plot1 <- ggplot(data = df[1:(eighty_percent + (0.25*eighty_percent)),], aes(x= reorder(PC, -Prop),y=Prop)) + 
    geom_bar(stat='identity') +
    geom_vline(aes(xintercept = eighty_percent), size = 2) +
    geom_label(aes(x = eighty_percent, y = 0.50,
                   label = '80%', vjust = -1, size = 8)) +
    geom_line(aes(y=cumulative,group=1 ), size = 2, color = 'red') +
    ggtitle('Variance explained by PCs') +
    xlab('Principal Component') +
    ylab('Proportion of variance explained (%)') +
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1), legend.position = 'none') +
    scale_x_discrete(breaks = df$PC[seq(1, length(df$PC), by = scale)])
  return(plot1)
}


#' Biplot of PCA
#'
#' @param pc PCA
#' @return biplot
make_biplot <- function(pc) {
  require(ggfortify)
  biplot1 <- autoplot(pc, loadings = TRUE) +
              xlab('PC1') + ylab('PC2') +
              theme(plot.title = element_text(size = 16, hjust = 0.5))
  return(biplot1)
}
  

#' Pairs plot
#'
#' @param pc PCA
#' @return pairsplot
make_pairsplot <- function(pc){
  plt <- ggplotify::as.ggplot(function() pairs(pc$rotation[,1:5], main = "Pairs plot"))
  return(plt)
}



#' Preprocessing of epigenetics metadata for heatmap
#'
#' @param unprocessed_metadata The epigenetics metadata
#' @return  Preprocessed metadata
preprocessing_epigenetics_metadata <- function(unprocessed_metadata) {
  unprocessed_metadata$Sex[which(is.na(unprocessed_metadata$Sex))] <- "M"
  df <- apply(unprocessed_metadata, 2, function(x) as.numeric(factor(x)))
  keep <- apply(df, 2, sd) > 0
  df <- df[, keep]
  keep <- !colSums(apply(df, 2, is.na)) == nrow(df)
  unprocessed_metadata <- scale(df[, keep])
  return(unprocessed_metadata[, -1]) # remove IDs
}


#' Perform correlation test with PCA values
#' 
#' @param pc The PCA
#' @param metadata Metadata that contain covariates for correlation testing
#' @param output_dir The directory where csv file is stored
#' @param name_file The name of the csv file
#' @return The correlation object
perform_PCA_correlation <- function(pc, metadata) {
  pc <- pc$x[,1:10]
  cxy <- cor(pc, metadata, use="pairwise.complete.obs")
  return(cxy)
}


#' Create heatmap of correlation of DNA methylation PCs and metadata.
#' 
#' @param cxy Correlation matrix.
#' @param output_dir Output directory for heatmap image file.
#' @return Heatmap plot object.
create_heatmap <- function(cor) {
  reshaped_correlation_matrix <- reshape2::melt(cor)
  plt <- ggplot(data = reshaped_correlation_matrix, aes(x=Var1, y=Var2, fill = value))+
          geom_tile() +
          ggtitle('Heatmap correlation matrix PCs and confounding factors') +
          xlab('Principal Components') +
          ylab('Confounding factors') +
          theme(plot.title = element_text(size = 16, hjust = 0.5))
  return(plt)
}


#' Subset numerical columns
#'
#' @param data input data
#' @return subset of numerical columns
subset_num_cols <- function(data) {
  nums <- unlist(lapply(data, is.numeric))
  subset <- data[,nums]
  return(subset)
}


#' Subset on matching Subjects using IDs file # Pre-processing for PCA correlation
#'
#' @param file1
#' @param file2
#' @param reference_file IDs file
#' @return subset of matching subjects
find_matching_IDs <- function(file1, file2, IDs) {

  large_file <- subset_num_cols(file1)
  values_file <- subset_num_cols(file2)
  
  IDs <- IDs[!is.na(IDs$XOmicsPhenoID) & !is.na(IDs$XOmicsmetaboID),]
  IDs <- subset(IDs, subset = IDs$XOmicsmetaboID %in% rownames(values_file))
  IDs <- subset(IDs, subset = IDs$XOmicsPhenoID %in% rownames(large_file))

  subset <- subset(large_file, subset = rownames(large_file) %in% IDs$XOmicsPhenoID)
  subset2 <- subset(values_file, subset = rownames(values_file) %in% IDs$XOmicsmetaboID)	
  
  return(list(subset, subset2))
}


#' Extract survey elements from phenotypes data
#'
#' @param phenotypes file
#' @return df with only survey elements
extract_survey_elements <- function(phenotypes) {
  survey_element_cols <- grep("q[0-9]+[a-z]+[0-9]+", colnames(phenotypes), value = TRUE)
  survey_colnames <- colnames(phenotypes)[which( !colnames(phenotypes) %in% survey_element_cols)]
  return(phenotypes[,survey_colnames])
}




# Targets

# Preprocessing of metabolomics data for PCA
metabolomics_preprocessed <- filter_samples_missing_values(metabolomics_values_normalized)


# Impute missing values
epigenomics_imputed <- data_imputation_median(epigenetics_preprocessed)


# Find matching samples between metabolomics and phenotypes data
matching_IDs_metabolomics_phenotypes <- find_matching_IDs(phenotypes_set2, metabolomics_preprocessed, IDs)


# Perform PCA
pca_list <- lapply(list(matching_IDs_metabolomics_phenotypes[[2]], t(epigenomics_imputed)), perform_PCA)


# Generate plots
screeplots <- lapply(pca_list, make_screeplot)
biplots <- lapply(pca_list, make_biplot)
pairsplots <- lapply(pca_list, make_pairsplot)


# Pre-processing epigenomics metada 
epigenetics_meta_preprocessed <- preprocessing_epigenetics_metadata(epigenetics_meta)

# Pre-processing of phenotypes data (metabolomics metadata)
phenotypes_survey_data <- extract_survey_elements(matching_IDs_metabolomics_phenotypes[[1]])

# Correlation of PCs with metadata
mtblmcs_pca_cor <- perform_PCA_correlation(pca_list[[1]], phenotypes_survey_data) 
epi_pca_cor <- perform_PCA_correlation(pca_list[[2]], epigenetics_meta_preprocessed)

# Create heatmaps
mtblmcs_heatmap <- create_heatmap(mtblmcs_pca_cor)
epi_heatmap <- create_heatmap(epi_pca_cor)

# Save all images to PDF
p <- list(screeplots,
          biplots,
          pairsplots,
          mtblmcs_heatmap,
          epi_heatmap)

pdf(output_dir)
p
dev.off()
