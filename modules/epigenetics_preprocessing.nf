#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir





////////////////////////////////////////////////////
/* --       Epigenomics preprocessing          -- */
////////////////////////////////////////////////////+

process EPIGENOMICS_ANNOTATION {

	label 'full_resources'	

	input:
	path epigenomics_MOESM1
	path epigenomics_MOESM4
	path epigenomics_MOESM5
	path MethylationEPIC
	path anno_epic

	output:
	path 'epigenomics_annotation_out.RData'

	"""
	Rscript $project_dir/Scripts/epigenomics_annotation.R ${epigenomics_MOESM1} ${epigenomics_MOESM4} ${epigenomics_MOESM5} ${MethylationEPIC} ${anno_epic}  'epigenomics_annotation_out.RData'
	"""
}



process EPIGENOMICS_FILTERING {
	
	label 'full_resources'
	
	input:
	path epigenomics_values
	path annotation_out

	output:
	path '*_filtered.csv'

	"""
	Rscript $project_dir/Scripts/epigenomics_filtering.R ${epigenomics_values} ${annotation_out}  '${epigenomics_values}_filtered.csv'
	"""
}


process EPIGENOMICS_IMPUTATION { 

	label 'full_resources'

	input:
	path epigenomics_values 

	output:
	path 'epigenomics_imputed_*.csv'
	
	"""
	Rscript $project_dir/Scripts/epigenomics_imputation.R ${epigenomics_values} 'epigenomics_imputed_${epigenomics_values}.csv'
	"""
}


process EPIGENOMICS_COVARIATES { 

	label 'full_resources'

	input:
	path epigenomics_values 
	path epigenomics_meta

	output:
	path 'epigenomics_corrected.csv'
	
	"""
	Rscript $project_dir/Scripts/CovariateCorrection.R ${epigenomics_values}  ${epigenomics_meta} 'epigenomics_corrected.csv'
	"""
}

process SUBSET_SD {

	label 'full_resources'

	input:
	path epigenomics_values
	val feature_cutoff

	output:
	path 'epigenomics_values_subset.csv'

	"""
	Rscript $project_dir/Scripts/sort_cols_sd.R ${epigenomics_values} ${feature_cutoff} epigenomics_values_subset.csv
	"""

}


process EPIGENOMICS_SCALING { 

	label 'full_resources'

	input:
	path epigenomics_values 

	output:
	path 'epigenomics_scaled.csv'
	
	"""
	Rscript $project_dir/Scripts/epigenomics_scaling.R ${epigenomics_values} 'epigenomics_scaled.csv'
	"""
}
