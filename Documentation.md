
# Input files

| Filename | Description  | Columns | Rows | File extension |
|---------|-------------|-------|--------| ------- |
| ids | Sample IDs for each data type  | data type  | sample  | `.csv`  |
| maf_files  |  [Metabolite Assignment File](https://www.ebi.ac.uk/metabolights/guides/MAF/Title)  | feature metadata/sample | metabolite | `.tsv` |
| behavioral_data| CBCL data | CBCL question | sample | `.csv` |
| epigenomics_meta | epigenomics sample metadata | sample metadata | sample  | `.csv` |
| epigenomics_values |  epigenomics data (bèta values)  | sample | CpG site | `.csv` |
| mtblmcs_values |  metabolite intensities | metabolite | sample  | `.csv` |
| phenotype_covariates | phenotypic data (used as covariates) | covariate | sample | `.csv`  |

# Workflow parameters

| Parameter | Description  | Default value | Possible values |
|---------|----------------| ------- | ----------|
| `missingness_cutoff` | Cutoff used to filter samples/feature on missing values   | 0.15 | 0-1|
| `mtblmcs_normalization_type` | Metabolites intensities normalization based on measured creatinine levels or specific gravity  | `cr` |  `cr`, `sg` | 
| `cbcl_imputation_method` | Imputation method for CBCL data. Method can be set to either `RF` (imputation by random forests using `missForest`) or `MCA` (imputation using MCA via `missMDA`). | `RF`| `RF`, `MCA` |
| `scale_epigenomics` | Epigenomics features are scaled by mean centering | false | false, true |
| `feature_subset_cutoff` | Cutoff used to subset epigenomics on percentage of features wiht highest variance | 10 | 0-100 |
| `seed` | Random seed for training the MOFA model | 100 | integer | 
| `convergence_mode` |   Argument to convergence mode while training the MOFA model | `slow`  | (`fast`, `medium`, `slow`) |
| `output`| Directory where results are stored | | |



# Workflow components

## Metabolomics pre-processing

| Process | Parameters (value options) | Description | Input | Output |
|---------|----------------------------|-------------|-------|--------|
| ` METABOLOMICS_FILTERING` | `missingness_cutoff` (Default = 0.15) | Rows and columns are filtered on missingness (in the respective order). A treshold of 15% is used.  |  `amines_MAF.tsv`, `OA_MAF.tsv`, `steroids_MAF.tsv`  | `metabolomics_filter.html`,  `*_filtered.csv` |
| `METABOLOMICS_NORMALIZATION`  | `mtblmcs_normalization_type` (`cr`, `sg`) | Performs metabolites intensities normalization based on measured creatinine levels ór specific gravity |  `*_filtered.csv`, mtblmcs_values  |  `*_normalized.csv` |
| `METABOLOMICS_SCALING`  |  | Scales metatobolite intensities using pareto scaling (scaled by the square-root of their column-wise standard deviations). Reference: https://doi.org/10.1186/1471-2164-7-142 |  `*_normalized.csv` |  `*_scaled.csv` |
| `CONCATENATE_METABOLOMICS` | | Concatenates the different metabolomics platforms, that were pre-processed separately |  `amines_scaled.csv`, `oa_scaled.csv`, `steroids_scaled.csv` | `metabolomics_preprocessed.csv` |

## Epigenomics pre-processing

| Process | Parameters (value options) | Description | Input | Output |
|---------|----------------------------|-------------|-------|--------|
| `EPIGENOMICS_ANNOTATION`  | | Annotate CpG sites using the R package 'meffil'. CpG IDs of only EPIC Probes and CpG IDs of both 450k and EPIC are identified. Secondly, CpG sites are annotated using the standard Illumina annotation file (methylationEPIC)   | Epigenomics MOESM1, MOESM4, MOESM5, methylationEPIC, epic_annotation | `epigenomics_annotation_out.RData` |
| `EPIGENOMICS_FILTERING`   | | Probes were excluded from all samples if they mapped to multiple locations in the genome, if they overlapped with a Single Nucleotide Polymorphism (SNP) or Insertion/Deletion (INDEL), or if they had a success rate < 0.95 across samples. Annotations of ambiguous mapping probes (based on an overlap of at least 47 bases per probe) and probes where genetic variants (SNPs or INDELS) with a minor allele frequency > 0.01 in Europeans overlap with the targeted CpG or single base extension site (SBE) were obtained from Pidsley et al (2016) (https://doi.org/10.1186/s13059-016-1066-1) | Epigenomics values file, `epigenomics_annotation_out.RData` | `epigenomics_values_filtered.csv` |
| `EPIGENOMICS_IMPUTATION`  |  |Missing values in the epigenomics beta values data were imputed using the median of the data. | `epigenomics_values_filtered.csv` | `epigenomics_imputed.csv` |
| `EPIGENOMICS_COVARIATES`  | | Correction of covariates is performed by regressing out the effect of sample plate number, sample row number and various cell type concentrations | `epigenomics_imputed.csv` | `epigenomics_corrected.csv` |
| `SUBSET_SD` | `feature_cutoff` (Default=10%) | The epigenomics values data is subsetted on the 10% features (CpG sites), showing the highest standard deviation. | `epigenomics_corrected.csv`  | `epigenomics_values_subset.csv` |
| `EPIGENOMICS_SCALING`  | | All values were scaled by mean centering. | `epigenomics_values_subset.csv` | `epigenomics_scaled.csv` |


## CBCL data pre-procssing and MCA

| Process | Parameters (value options) | Description | Input | Output |
|---------|----------------------------|-------------|-------|--------|
| `CBCL_FILTER_IMPUTE_MCA` | `cbcl_imputation_method` (`RF`, `MCA`) | Process filters selected phenotype data based on missingness, imputes missing values, and performs Multiple Correspondance Analysis (MCA) (DOI: 10.1285/i20705948v10n2p432). Imputation method can be chosen by setting `cbcl_imputation_method` to either `RF` (imputation by random forests using `missForest`) or `MCA` (imputation using MCA via `missMDA`). Individuals' coordinates for new MCA dimensions (principal coordinates) are stored in `cbcl_mca_coord.csv` with individuals in rows and dimensions in columns. | `cbcl_data.csv` | `CBCL_filter_impute_MCA.html`, `cbcl_filtered.csv`, `cbcl_imputed.csv`, `cbcl_mca_coord.csv`, `cbcl_mca.RData` |

## Omics data sample identifier mapping
| Process | Description | Input | Output |
|---------|-------------|-------|--------|
| `MAP_IDS`  | The epigenomics and metabolomics preprocesssed data do not contain the same samples. Here, the respective dataframes are subsetted on the samples that can be found in both preprocessed datasets. These mapped dataframes can be used for the integration methods (SNF, MOFA) | `epigenomics_values_subset.csv`,  `metabolomics_preprocessed.csv`, IDs file | `epigenomics_values_mapped.csv`, `metabolomics_values_mapped.csv`, `duplicates_discarded.csv` |

## Analysis
| Process | Parameters (value options) | Description | Input | Output |
|---------|----------------------------|-------------|-------|--------|
| `HEATMAP_MISSINGNESS`   |  | Heatmaps showing missing values in all the datasets -  as uploaded to the DRE - are generated.  | Metabolomics values file, Metabolomics metadata file, Epigenomics values file, Epigenomics metadata file, Phenotypes SPSS file | `heatmap_missingness.pdf` |
| `PCA`   |  | Principal Component Analysis was performed on the metabolomics and epigenomics data, using the prcomp() function in R. Screeplots, biplots and pairsplots were generated of the different PCAs. Correlation of the first 10 PCs with metadata (in Phenotypes SPSS file and epigenomics metdata) was vizualized in heatmaps.   | `metabolomics_corrected.csv`,  `epigenomics_scaled.csv`, Epigenomics metadata file,  `cbcl_imputed.csv`, IDs file | `pca.pdf` |
| `SNF`   |  | Similarity Network Fusion (Wang, B., Mezlini, A., Demir, F. et al. Similarity network fusion for aggregating data types on a genomic scale. Nat Methods 11, 333–337 (2014) (https://doi.org/10.1038/nmeth.2810) is performed on the epigenomics and metabolomics data.  | `epigenomics_values_mapped.csv`, `metabolomics_values_mapped.csv` | `snf.pdf`, `snf.npy` |
| `SNF_ANALYSIS`   |  | The sample correlation matrix produced by SNF is analyzed using spectral clustering. Then, clusters are analyzed for possible differences in phenotype covariates and behavioral dimensions (MCA)  | `snf.csv`, Phenotype covariates, `metabolomics_values_mapped.csv`  | `snf_analysis_out.ipynb` |
| `MOFA`  | `seed` (integer) and `convergence_mode` (`fast`, `medium`, `slow`) | Argelaguet, R. et al. Multi-Omics Factor Analysis—a framework for unsupervised integration of multi-omics data sets. Mol Syst Biol. 14:e812 (2018) (https://doi.org/10.15252/msb.20178124) is performed on the epigenomics and metabolomics data. | `epigenomics_values_mapped.csv`, `metabolomics_values_mapped.csv` | `MOFAmodel.hdf5`  |
