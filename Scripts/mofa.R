#!/usr/bin/Rscript

################################################################################
## Title:         mofa.R
## Description:   Perform multi-omics factor analysis (MOFA) on a list of omics data matrices
## Author:        Purva Kulkarni
## Date created:  2022-15-06
## Email:         purva.kulkarni@radboudumc.nl
################################################################################
## Notes:
## 
################################################################################

# Load necessary libraries

library(MOFA2)
library(data.table)
library(ggplot2)
library(tidyverse)
library(ggpubr)


#Read in data

args = commandArgs(trailingOnly=TRUE)

number_of_arguments <- length(args)
number_of_omics <- number_of_arguments-3


#' Function to read in the omics data, transpose rows and columns
#'
#' @param omics_input path to omics data
#' @return data matrix
read_omics_data <- function(omics_input) {
	data <- read.csv(omics_input)
	#Transpose the data to have samples in columns and measured omics features in rows
	transposed_data <- as.matrix(setNames(data.frame(t(data[,-1])), data[,1]))
	return(transposed_data)
}

print(args)
print(args[1:number_of_omics])

# Argument for input data in the form of a list of matrices
omicsData <- lapply(args[1:number_of_omics], function(x) {read_omics_data(x)})
print(omicsData)
names(omicsData) <- c("epigenomics", "metabolomics") #TODO: Make these input parameters, and set default to 'omics1', 'omics2' etc.

convergence_mode = args[number_of_arguments-2] # Argument to convergence mode while training model :"fast", "medium", "slow"
print(convergence_mode)
seed = args[number_of_arguments-1] # Random seed for training the model
output_dir = args[number_of_arguments]

#' Function to perform MOFA
#' Input: A list containing multiple omics matrices
#' 
#' @param omicsData - List of omics data matrices with samples in columns and variables in rows
#' @convergence_mode -  Argument to convergence mode while training model :"fast", "medium", "slow"
#' @param seed - Random seed to train model
#' @return standard deviation 
#` @return output directory where the computed MOFA model file will be saved
performMOFA <- function(omicsData, convergence_mode, seed, output_dir) {
 
   # Create MOFA object and train the model
  MOFAobject <- create_mofa(omicsData)
  
  # Define data options
  # Important arguments:
  # scale_groups: if groups have different ranges/variances, it is good practice to scale each group to unit variance. Default is FALSE
  # scale_views: if views have different ranges/variances, it is good practice to scale each view to unit variance. Default is FALSE
  
  data_opts <- get_default_data_options(MOFAobject)
  
  # Define model options
  # Important arguments
  # num_factors: number of factors
  # likelihoods: likelihood per view (options are “gaussian”, “poisson”, “bernoulli”). Default is “gaussian”.
  # spikeslab_factors: use spike-slab sparsity prior in the factors? Default is  FALSE.
  # spikeslab_weights: use spike-slab sparsity prior in the weights? Default is  TRUE.
  # ard_factors: use ARD prior in the factors? Default is TRUE if using multiple groups.
  # ard_weights: use ARD prior in the weights? Default is TRUE if using multiple views.
  # Only change the default model options if familiar with the underlying mathematical model.
  
  model_opts <- get_default_model_options(MOFAobject)
  model_opts$num_factors <- 10
  
  # Define training options
  # Important arguments
  # maxiter: number of iterations. Default is 1000.
  # convergence_mode: “fast” (default), “medium”, “slow”. For exploration, the fast mode is sufficient. For a final model, consider using medium" or even “slow”, but hopefully results should not change much.
  # gpu_mode: use GPU mode? (needs cupy installed and a functional GPU).
  # verbose: verbose mode?
  # seed: random seed
  
  train_opts <- get_default_training_options(MOFAobject)
  train_opts$convergence_mode <- convergence_mode
  train_opts$seed <- seed
  
  # Prepare and train the MOFA object
  
  MOFAobject <- prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  
  # Train the model and save the computed model in a hdf5
  # outfile = file.path(output_dir,"MOFAmodel.hdf5")
  MOFAobject.trained <- run_mofa(MOFAobject, output_dir)
  
  return(MOFAobject.trained)
}

mofaobject_trained <- performMOFA(omicsData, convergence_mode, seed, output_dir) 
