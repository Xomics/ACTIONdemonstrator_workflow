#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir



process MOFA {

	publishDir "${params.output}/MOFA", mode: 'copy', overwrite: true

	input:
	path omicsDataList
	val seed
	val convergence_mode


	output:
	path 'MOFAmodel.hdf5'

	"""
	Rscript $project_dir/bin/mofa.R ${omicsDataList} ${convergence_mode} ${seed} MOFAmodel.hdf5
	"""

}

