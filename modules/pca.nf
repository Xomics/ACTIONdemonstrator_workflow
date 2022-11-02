#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir


////////////////////////////////////////////////////
/* --       Principal Component Analysis       -- */
////////////////////////////////////////////////////+


process PCA {

	label 'full_resources'

	publishDir "${params.output}/PCA", mode: 'copy', overwrite: true

	input: 
	path mtblmcs_values
	path epigenomics_values
	path epigenomics_meta
	path phenotypes_set2
	path ids
	
	output:
	path 'pca.pdf' 

	"""
	Rscript $project_dir/Scripts/pca.R ${mtblmcs_values} ${epigenomics_values} ${epigenomics_meta} ${phenotypes_set2} ${ids} pca.pdf
	"""
}
