# Dockerfile VERSION = v0.4
# docker login registry.cmbi.umcn.nl
# docker build -t registry.cmbi.umcn.nl/x-omics-action-dataset/action_nextflow/r-base-epigenomics-pre:$VERSION . 
# docker push registry.cmbi.umcn.nl/x-omics-action-dataset/action_nextflow/r-base-epigenomics-pre:$VERSION
# sudo docker pull registry.cmbi.umcn.nl/x-omics-action-dataset/action_nextflow/r-base-epigenomics-pre$VERSION
# sudo docker images # to get IMAGE_ID
# sudo docker save $IMAGE_ID -o r-base-epigenomics-pre.tar
# sudo singularity build r-base-epigenomics-pre.sif docker-archive://r-base-epigenomics-pre.tar
FROM r-base:4.1.2

# Installations needed for devtools package
RUN apt update
RUN apt install net-tools -y
RUN apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev -y

# Installations needed for R Cairo package
RUN apt-get install libcairo2-dev -y
RUN apt-get install libxt-dev -y

# Installations needed for R lme4 package
RUN apt install cmake -y
#RUN apt install libudunits2-dev -y

# Installation of BiocManager and biocLite
RUN R -e "install.packages('BiocManager', version = '3.12')" \
        && R -e "BiocManager::install(version = '3.14')" \ 
        && R -e "BiocManager::install('biocLite')" 

# R installations needed for devtools and devtools itself
RUN R -e "install.packages('installr')" \
        && R -e "installr::install.Rtools(check = TRUE, check_r_update = FALSE, GUI = FALSE)" \
        && R -e "install.packages(c('systemfonts', 'textshaping', 'ragg', 'pkgdown', 'Cairo', 'lme4'))" \
        && R -e "install.packages('devtools')" 


# Installation of meffil package from GitHub (v1.1.1) 
RUN R -e "devtools::install_github('perishky/meffil')"

# Installations for Nextflow metrics, 'ps' command
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*