#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir



process MAP_IDS {

	label 'full_resources'

	input:
	path epigenomics_values
	path mtblmcs_values
	path ids

	output:
	tuple val(1), path('epigenomics_values_mapped.csv')
	tuple val(1), path('metabolomics_values_mapped.csv')
	path 'duplicates_discarded.csv'

	"""
	python3 $project_dir/bin/map_IDs.py ${epigenomics_values} ${mtblmcs_values} ${ids} epigenomics_values_mapped.csv metabolomics_values_mapped.csv duplicates_discarded.csv
	"""

}



