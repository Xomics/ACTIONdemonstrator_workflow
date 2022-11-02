#!/usr/bin/env nextflow


nextflow.enable.dsl=2

project_dir = projectDir

////////////////////////////////////////////////////
/*    --               Functions               -- */
////////////////////////////////////////////////////+


def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run action.nf 
            --output dir/of/choice
            --mtblmcs_values Sythetic_data/synthetic_metabolomics.csv
            --maf_files Synthetic_data/
            --epigenomics_values Synthetic_data/synthetic_epigenomics.csv
            --epigenomics_meta Synthetic_data/synthetic_epigenomics_meta.csv
            --behavioral_data Synthetic_data/synthetic_cbcl_data.csv
            --phenotype_covariates Synthetic_data/synthetic_phenotype_covariates_data.csv
            --ids Synthetic_data/ACTIONdemonstrator_XOmics_IDS_synthetic.csv

       Optional arguments:
	     --container_dir                The directory where the required Singularity (.sif files) images are stored
         --mtblmcs_normalizaton_type    Normalization of metabolomics is by default done using creatinine (`cr`). Set to 'sg' for normalization by specific gravity.
         --cbcl_imputation_method		Behavioral data (CBCL) is imputed by default with random forests (`RF`). Set to 'MCA' for imputation with MCA.
         --covergence_mode				Argument to convergence mode while training MOFA model :"fast", "medium", "slow"
         --seed							Random seed to train MOFA model
         --feature_subset_cutoff		Percentage of features to select for analysis (for big dataframes). Features with highest standard deviation are selected.
		 --scale_epigenomics			Scale the epigenomics beta values by mean centering and dividing columns by standard deviation
        """
}



////////////////////////////////////////////////////
/* --            Input data files              -- */
////////////////////////////////////////////////////+


mtblmcs_values = Channel.fromPath("${params.mtblmcs_values}")  
maf_files = Channel.fromPath("${params.maf_files}/*_MAF.tsv")

phenotype_covariates = Channel.fromPath("${params.phenotype_covariates}") 
behavioral_data = Channel.fromPath("${params.behavioral_data}") 

epigenomics_values = Channel.fromPath("${params.epigenomics_values}")  
epigenomics_meta = Channel.fromPath("${params.epigenomics_meta}")   

epigenomics_MOESM1 = Channel.fromPath("EPIC_annotation/raw/13059_2016_1066_MOESM1_ESM.csv")
epigenomics_MOESM4 = Channel.fromPath("EPIC_annotation/raw/13059_2016_1066_MOESM4_ESM.csv")
epigenomics_MOESM5 = Channel.fromPath("EPIC_annotation/raw/13059_2016_1066_MOESM5_ESM.csv")
methylationEPIC =  Channel.fromPath("EPIC_annotation/raw/MethylationEPIC_v-1-0_B4.csv")
epic_annotation =  Channel.fromPath("EPIC_annotation/anno_epic_072017.RData")

ids = Channel.fromPath("${params.ids}") 



////////////////////////////////////////////////////
/* --                  Modules                 -- */
////////////////////////////////////////////////////+


include { METABOLOMICS_FILTERING; METABOLOMICS_NORMALIZATION; METABOLOMICS_SCALING } from './modules/metabolomics_preprocessing'
include { CONCATENATE_METABOLOMICS } from './modules/metabolomics_preprocessing'
include { EPIGENOMICS_ANNOTATION; EPIGENOMICS_FILTERING; EPIGENOMICS_IMPUTATION; SUBSET_SD; EPIGENOMICS_COVARIATES; EPIGENOMICS_SCALING } from './modules/epigenetics_preprocessing'
include { CBCL_FILTER_IMPUTE_MCA } from './modules/CBCL_MCA'
include { MAP_IDS } from './modules/map_IDs'
include { HEATMAP_MISSINGNESS } from './modules/heatmap_missingness'
include { PCA } from './modules/pca'
include { SNF; SNF_ANALYSIS } from './modules/snf'
include { MOFA } from './modules/mofa'



////////////////////////////////////////////////////
/* --                 Functions               -- */
////////////////////////////////////////////////////+


def group_mapped_omics(omics1, omics2) {

	omics_list = omics1
			.join(omics2)
	omics_list = omics_list.collect()
	return(omics_list.minus([1]))

}



////////////////////////////////////////////////////
/* --                 Workflow                 -- */
////////////////////////////////////////////////////+


workflow {

	// Show help message
	if (params.help) {
    	helpMessage()
    	exit 0
	}


	////////////////
	// SUBWORKFLOW 1: Metabolomics data pre-processing
	////////////////
	METABOLOMICS_FILTERING(maf_files, params.missingness_cutoff) 
	METABOLOMICS_NORMALIZATION(METABOLOMICS_FILTERING.out[1].flatten(), mtblmcs_values.toList(), params.mtblmcs_normalization_type) 
	METABOLOMICS_SCALING(METABOLOMICS_NORMALIZATION.out.flatten())
	CONCATENATE_METABOLOMICS(METABOLOMICS_SCALING.out.collect())


	////////////////
	// SUBWORKFLOW 2: Epigenomics data pre-processing
	////////////////
	EPIGENOMICS_ANNOTATION(epigenomics_MOESM1, epigenomics_MOESM4, epigenomics_MOESM5, methylationEPIC, epic_annotation)
	EPIGENOMICS_FILTERING(epigenomics_values, EPIGENOMICS_ANNOTATION.out)
	EPIGENOMICS_IMPUTATION(EPIGENOMICS_FILTERING.out)
	EPIGENOMICS_COVARIATES(EPIGENOMICS_IMPUTATION.out, epigenomics_meta)
	SUBSET_SD(EPIGENOMICS_COVARIATES.out, params.feature_subset_cutoff)
	EPIGENOMICS_SCALING(SUBSET_SD.out)


	////////////////
	// SUBWORKFLOW 3: Behavior data imputation and MCA
	////////////////
	CBCL_FILTER_IMPUTE_MCA(behavioral_data)


	////////////////
	// SUBWORKFLOW 4: Map omics data on sample IDs
	////////////////
	MAP_IDS(SUBSET_SD.out, CONCATENATE_METABOLOMICS.out, ids)


	////////////////
	// SUBWORKFLOW 5: Splitting csv file in chunks for lower memory usage
	////////////////
	//chunksize = Channel.value(10000)
	// EPIGENOMICS_SPLIT_CHUNKS(epigenomics_values_head, chunksize)
	// EPIGENOMICS_FILTERING(EPIGENOMICS_SPLIT_CHUNKS.out.flatten(), epigenomics_annotation.out.toList())
	// CONCATENATE_EPI_CHUNKS(EPIGENOMICS_FILTERING.out.toList())


	////////////////
	// SUBWORKFLOW 6: Create heatmap of missingness (unprocessed data)
	////////////////
	HEATMAP_MISSINGNESS(mtblmcs_values, epigenomics_values, epigenomics_meta, behavioral_data, phenotype_covariates) 


	////////////////
	// SUBWORKFLOW 7: Principal Component Analysis
	////////////////
	PCA(CONCATENATE_METABOLOMICS.out, EPIGENOMICS_SCALING.out, epigenomics_meta, CBCL_FILTER_IMPUTE_MCA.out[2], ids)


	////////////////
	// SUBWORKFLOW 8: Multi-omics Factor Analysis
	////////////////
	omics_list = group_mapped_omics(MAP_IDS.out[0], MAP_IDS.out[1])
	MOFA(omics_list, params.seed, params.convergence_mode)


	////////////////
	// SUBWORKFLOW 9: Similarity Network Fusion
	////////////////
	SNF(omics_list)
	SNF_ANALYSIS(SNF.out[1], phenotype_covariates, MAP_IDS.out[3], CBCL_FILTER_IMPUTE_MCA.out[3], params.output)
}
