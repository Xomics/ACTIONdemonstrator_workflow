# Dockerfile VERSION = v0.4
# docker login registry.cmbi.umcn.nl
# docker build -t registry.cmbi.umcn.nl/x-omics-action-dataset/action_nextflow/r-base-analysis:$VERSION . 
# docker push registry.cmbi.umcn.nl/x-omics-action-dataset/action_nextflow/r-base-analysis:$VERSION
# sudo docker pull registry.cmbi.umcn.nl/x-omics-action-dataset/action_nextflow/r-base-analysis$VERSION
# sudo docker images # to get IMAGE_ID
# sudo docker save $IMAGE_ID -o r-base-analysis.tar
# sudo singularity build r-base-analysis.sif docker-archive://r-base-analysis.tar
FROM r-base:4.1.2

# Install cmake and GSL
RUN  apt-get update && apt-get install -y cmake \
    &&  apt-get install -y libgsl-dev

# Install R packages
RUN R -e "install.packages(c('nlme', 'foreach', 'doMC', 'factoextra',  'ggplot2' ), dependencies = TRUE, repos = 'http://cran.rstudio.com/')" 
RUN R -e "install.packages(c('foreign', 'readxl', 'visdat', 'ggplotify'), dependencies = TRUE, repos = 'http://cran.rstudio.com/')" 
RUN R -e "install.packages(c('reshape2', 'gridExtra', 'ggfortify', 'IMIFA'), dependencies = TRUE, repos = 'http://cran.rstudio.com/')"

# Installations for Nextflow metrics, 'ps' command
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*