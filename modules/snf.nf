#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir





process SNF {

	label 'full_resources'

	publishDir "${params.output}/SNF", mode: 'copy', overwrite: true
	
	input: 
	path omicsData

	output:
	path 'snf.csv'

	"""
	python3 $project_dir/bin/perform_snf.py ${omicsData}  snf.csv 
	"""
}


process SNF_ANALYSIS {

	publishDir "${params.output}/SNF", mode: 'copy', overwrite: true

	input:
	path snf_array
	path phenotype_covariates
	path mca_dims
	path output_dir_plots

	output:
	path 'snf_analysis_out.ipynb'
	path 'snf_analysis_out.csv'	
	
	"""
	python -m ipykernel install --name base --user
	papermill $project_dir/bin/snf_analysis.ipynb  snf_analysis_out.ipynb -p snf_matrix_path ${snf_array} -p phenotypes_covariates_path ${phenotype_covariates} -p mca_dims_path ${mca_dims} -p output_dir_plots ${output_dir_plots} -p outdir_df 'snf_analysis_out.csv' -k base 
	"""

}

process SNF_GEE {
	
	publishDir "${params.output}/SNF", mode: 'copy', overwrite: true

	input:
	path snf_analysis_out
	
	output:
	path 'snf_gee.html'

	"""
	cp -L $project_dir/bin/snf_gee_analysis.Rmd snf_gee_analysis.Rmd
	Rscript -e "rmarkdown::render('snf_gee_analysis.Rmd', output_format = 'html_document', output_file = 'snf_gee.html', params = list(input_file_snf_pheno = '${snf_analysis_out}' ))"
	"""



}
