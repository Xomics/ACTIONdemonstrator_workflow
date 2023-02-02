#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir




////////////////////////////////////////////////////
/* --       Phenotypes imputation and MCA       -- */
////////////////////////////////////////////////////+


process CBCL_FILTER_IMPUTE_MCA { 
	
	input:
	path pheno
	path cbcl_labels_path

	output:
	path 'CBCL_filter_impute_MCA.html' 
	path 'cbcl_filtered.csv' 
	path 'cbcl_imputed.csv' 
	path 'cbcl_mca_coord.csv' 
	path 'cbcl_mca.RData' 
	
	"""
	cp -L $project_dir/bin/CBCL_filter_impute_MCA.Rmd CBCL_filter_impute_MCA.Rmd
	Rscript -e "rmarkdown::render('CBCL_filter_impute_MCA.Rmd', output_format = 'html_document', output_file = 'CBCL_filter_impute_MCA.html', params = list(cbcl_infile = '${pheno}', cbcl_labels = '${cbcl_labels_path}', cbcl_filtered_outfile = 'cbcl_filtered.csv', cbcl_imputed_outfile = 'cbcl_imputed.csv', cbcl_mca_coord_outfile = 'cbcl_mca_coord.csv', cbcl_mca_outfile = 'cbcl_mca.RData', imputation_method = '${params.cbcl_imputation_method}'))"	
	"""
}



