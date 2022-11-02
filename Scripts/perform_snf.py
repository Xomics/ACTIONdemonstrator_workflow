import numpy as np
import sys
import os
import snf
import pandas
from snf import cv
from sklearn.cluster import spectral_clustering
from sklearn.metrics import v_measure_score


number_of_arguments = len(sys.argv)

output_pdf = sys.argv[number_of_arguments-2]
output_npy = sys.argv[number_of_arguments-1]

def read_csv_file(file_path):
	df = pandas.read_csv(file_path, index_col=0)
	#Remove NA values
	df = df.dropna()
	return(df)

def find_common_IDs(df_list):
	#find common IDS
	index_list = []
	for i in df_list:
		index_list.append(i.index)
	common_IDs = set.intersection(*map(set,index_list))

	#subset every df on these IDS
	new_df_list = []
	for i in df_list:
		df_with_common_ids = i[i.index.isin(common_IDs)]	
		array = df_with_common_ids.to_numpy() #convert to Numpy Array
		new_df_list.append(array)
	
	return(new_df_list)


# Read in files (independent of number of input files)
df_list = []
for n in range(1,number_of_arguments-2):
	df = read_csv_file(sys.argv[n])
	df_list.append(df)


# Subset on common IDs and convert to Numpy Array
data_combined = find_common_IDs(df_list)

# Cross validation
#grid_zaff, grid_labels = cv.snf_gridsearch(data_combined)
#mu_cv, k_cv = cv.get_optimal_params(grid_zaff, grid_lables)


# Calculate affinity networks for single omics
affinity_networks = snf.make_affinity(data_combined, metric = 'euclidean', K=20, mu=0.5)


# Remove infinite values, set to 0 
for i in affinity_networks:
	i[~np.isfinite(i)] = 0


# SNF
fused_network = snf.snf(affinity_networks)


# Safe fused network in binary format
np.savetxt(output_npy, fused_network, delimiter=",")


# Save cross validation results to PDF
with open (output_pdf, 'w') as f:
	f.write('Cross validaton results for K and mu')
#	f.write('K = {:.2f}'.format(k_cv))
#	f.write('mu = {:.2f}'.format(mu_cv))
