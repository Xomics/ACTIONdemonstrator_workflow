#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir



////////////////////////////////////////////////////
/* --            Heatmap missingness            -- */
////////////////////////////////////////////////////+


process HEATMAP_MISSINGNESS {

	label 'r_base_analysis_small_tasks'

	publishDir "${params.output}/Missing_data", mode: 'move', overwrite: true

	input: 
	path mtblmcs_values
	path epigenomics_values
	path epigenomics_meta
	path behavioral_data
	path phenotype_covariates

	output:
	path 'heatmap_missingness.pdf' 

	"""
	Rscript $project_dir/bin/heatmap_missingness.R ${mtblmcs_values} ${epigenomics_values} ${epigenomics_meta} ${behavioral_data} ${phenotype_covariates} heatmap_missingness.pdf
	"""
}

