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

process MOFA_ANALYSIS {

	publishDir "${params.output}/MOFA", mode: 'copy', overwrite: true

	input:
	path mofa_model
	path pheno_covariates
	path cbcl_mca_dimensions

	output:
	path 'mofa_analysis.html'

	"""
	cp -L $project_dir/bin/MOFA_downstream_analysis_report.Rmd MOFA_downstream_analysis_report.Rmd
	Rscript -e "rmarkdown::render('MOFA_downstream_analysis_report.Rmd', output_format = 'html_document', output_file = 'mofa_analysis.html',  params = list(modelFile = '${mofa_model}', metaData = '${pheno_covariates}', mcaCoord = '${cbcl_mca_dimensions}'  ))"
	"""

}

