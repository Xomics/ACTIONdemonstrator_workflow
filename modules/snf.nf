#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir





process SNF {

	label 'full_resources'

	publishDir "${params.output}/SNF", mode: 'copy', overwrite: true
	
	input: 
	path omicsData

	output:
	path 'snf.pdf' 
	path 'snf.csv'

	"""
	python3 $project_dir/bin/perform_snf.py ${omicsData} snf.pdf snf.csv 
	"""
}


process SNF_ANALYSIS {

	publishDir "${params.output}/SNF", mode: 'copy', overwrite: true

	input:
	path snf_array
	path phenotype_covariates
	path metabolomics
	path mca_dims
	path output_dir_plots

	output:
	path 'snf_analysis_out.ipynb'

	"""
	python -m ipykernel install --name base --user
	papermill $project_dir/bin/snf_analysis.ipynb  snf_analysis_out.ipynb -p snf_matrix_path ${snf_array} -p phenotypes_covariates_path ${phenotype_covariates} -p metabolomics_path ${metabolomics} -p mca_dims_path ${mca_dims} -p output_dir_plots ${output_dir_plots} -k base
	"""

}
