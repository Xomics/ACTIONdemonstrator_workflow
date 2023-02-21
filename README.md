# X-omics ACTIONdemonstrator multi-omics analysis workflow

Here, the ACTIONdemonstrator workflow is made using Nextflow workflow manager. With this workflow, we would like to demonstrate how to create a FAIR workflow. By making this a containerzied/modularized workflow and compatible with WorkflowHub, we consider this a FAIR workflow. 


# Pipeline summary

The NTR-ACTION analysis workflow consists of a set of different pre-processing and analysis steps:

![alt text](flowchart.png)
[More complete documentation on the workflow can be found here](Documentation.md)


# WorkflowHub and RO-Crate

The workflow can be found on the  [WorkflowHub page](https://workflowhub.eu/workflows/402) page as well. Different versions of the workflow will be published here. 
The webpage provides an option to download the workflow as Research-Object Crate (RO-Crate). This RO-Crate contains all required workflow elements - including metadata describibing these elements - packed into one `.zip` file. The Jupyter Notebook [ `generate_ro-crate.ipynb`](generate_ro-crate.ipynb) demonstrates how the RO-Crate was made. 


# Requirements

Nextflow and Singularity are required to be installed. To install them on CentOS Linux 7 enable [EPEL](https://docs.fedoraproject.org/en-US/epel/) and run

```
sudo yum -y install conda.noarch
sudo yum -y install singularity.x86_64
```

## Installing Nextflow
```
conda create --name nf
conda activate nf
conda install -c bioconda nextflow:22.04.0
```

See also https://anaconda.org/bioconda/nextflow and https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html.


# Usage

## Obtain workflow

The workflow RO-Crate can be downloaded from the [WorkflowHub page](https://workflowhub.eu/workflows/402). The RO-Crate contains all necessary files to run the workflow (where this repository does not because of storage limits).

Alternatively, you can clone this git repository (`git clone https://github.com/Xomics/ACTIONdemonstrator_workflow.git`). Then download the missing Infinium MethylationEPIC product file [here](https://emea.support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html), for which a placeholder is located in `EPIC_annotation/raw/`.


## Run main workflow:

For detailed instructions on how to run the workflow:
```
nextflow run action.nf --help
```

The typical command to run the workflow is:
```
nextflow run action.nf 
            --output dir/of/choice
            --mtblmcs_values Sythetic_data/synthetic_metabolomics.csv
            --maf_files Synthetic_data/
            --epigenomics_values Synthetic_data/synthetic_epigenomics.csv
            --epigenomics_meta Synthetic_data/synthetic_epigenomics_meta.csv
            --behavioral_data Synthetic_data/synthetic_cbcl_data.csv
            --phenotype_covariates Synthetic_data/synthetic_phenotype_covariates_data.csv
            --cbcl_labels Synthetic_data/cbcl_labels.csv
            --ids Synthetic_data/ACTIONdemonstrator_XOmics_IDS_synthetic.csv
```
To run the workflow **in the [Digital Research Environment](https://mydre.org/)**:
- Specify config file: `-c dre.config`
- Specify container files (.sif) dir: `--container_dir /dir/to/containers/`



# Synthetic Data
Synthetic data is provided as example data in `Synthetic_data/`. This directory contains the necessary (omics) data files to run the workflow + the ISA-Tab files describing the experimental metadata. These ISA files were generated with this [repository](https://github.com/Xomics/ISA-ACTION-Template).

# Software containers
Nextflow automatically pulls all necessary Docker containers from Dockerhub, when using the default `nextflow.config` file. Docker (and Singularity) containers can be built locally as well, see more detailed information [in this repository](https://github.com/Xomics/Docker_containers.git). The `dre.config` file demonstrates how to allocate local Singularity (.sif) files to Nextflow processes. 


# Authors

Radboud University Medical Center, Nijmegen, Netherlands:
- Anna Niehues 
- Casper de Visser
- Purva Kulkarni 
- Alain J. van Gool 
- Peter A. C. 't Hoen

Vrije Universiteit Amsterdam, Amsterdam, Netherlands:
- Fiona A. Hagenbeek
- René Pool 
- Dorret I. Boomsma 
- Jenny van Dongen 

Leiden University, Leiden, Netherlands:
- Naama Karu 
- Alida S. D. Kindt 