#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir



process METABOLOMICS_FILTERING { 

	input:
	path MAF
	val missingness_cutoff

	output:
	path 'metabolomics_filter.html' 
	path '*_filtered.csv'
	
	"""
	cp -L $project_dir/Scripts/metabolomics_filter.Rmd metabolomics_filter.Rmd
	Rscript -e "rmarkdown::render('metabolomics_filter.Rmd', output_format = 'html_document', output_file = 'metabolomics_filter.html',  params = list(MAF_path = '${MAF}', output_dir_csv = '${MAF}_filtered.csv', missing_value_cutoff = '${missingness_cutoff}'  ))"
	"""
}

process METABOLOMICS_NORMALIZATION { 	

	label 'r_base_analysis_small_tasks'
	
	input:
	path mtblmcs_filtered 
	path mtblmcs_values
	val normalization_type
	
	output:	
	path '*_normalized.csv'
	
	script:
	"""
	Rscript $project_dir/Scripts/metabolomics_normalization.R ${mtblmcs_filtered} ${mtblmcs_values} ${mtblmcs_filtered}_normalized.csv ${normalization_type}
	"""
}

process METABOLOMICS_SCALING { 	

	label 'r_base_analysis_small_tasks'
	
	input: 
	path mtblmcs_values
	
	output:	
	path '*_scaled.csv'
	
	script:
	"""
	Rscript $project_dir/Scripts/metabolomics_scaling.R  ${mtblmcs_values} ${mtblmcs_values}_scaled.csv 
	"""
}

process CONCATENATE_METABOLOMICS {

	label 'r_base_analysis_small_tasks'

 	input:
 	path mtblmcs_files
  
	output:
	path 'metabolomics_preprocessed.csv'
  
	script:
	"""
	Rscript $project_dir/Scripts/concatenate_MAF.R ${mtblmcs_files} metabolomics_preprocessed.csv
	"""

}



