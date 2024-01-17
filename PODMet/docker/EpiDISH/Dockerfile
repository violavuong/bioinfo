# Build EpiDISH
## Base image
FROM bioconductor/bioconductor_docker:3.17

## Data/time of build
ARG buildDate

## Requirements needed to run EpiDISH R script
ADD EpiDISH.R /source/

######### METADATA #########
# Image information
LABEL about.image="Bioconductor base image with required libraries for EpiDISH"
LABEL image.authors="Viola Meixian Vuong"
LABEL image.authors.email="ViolaMeixian.Vuong@UGent.be"
## Update when new R and libraries version are available
LABEL image.maintainer="Maisa R Ferro dos Santos"
LABEL image.maintainer.email="MaisaRenata.FerrodosSantos@ugent.be"
LABEL base.image="bioconductor/bioconductor_docker:3.17"
LABEL image.version="0.1"
## Date/time of build
LABEL buildDate=${buildDate}
# Software information
LABEL software="EpiDISH for DNA methylation deconvolution in EWAS studies"
LABEL software.version="2.16.0"
LABEL about.article="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1511-5"
LABEL about.documentation="https://www.bioconductor.org/packages/release/bioc/html/EpiDISH.html, https://github.com/sjczheng/EpiDISH"
LABEL about.tags="DNA methylation, EWAS, Deconvolution"
# Docker command examples for running this image
LABEL docker.run01.help="docker run --rm epidish:0.1 /source/EpiDISH.R -h"
LABEL docker.run01.deconvolution="docker run --rm -v path-to-data-directory:/mnt -it epidish:0.1 /source/EpiDISH.R /mnt/path-to-ref /mnt/path-to-samples modality /mnt/path-to-outputs"
###########################################################

## Check for updates
RUN apt-get update -y && apt-get install -y\
 libcurl4-gnutls-dev \
 libssl-dev \
 libssh2-1-dev \
 libxml2-dev \
 libhdf5-dev \
 libgmp-dev \
 libpng-dev \
 libgsl-dev \
 libxt-dev \
 libcairo2-dev \
 libtiff-dev \
 fftw-dev \
 procps

## Install EpiDISH dependencies
RUN R -e 'BiocManager::install(c("e1071","quadprog","matrixStats","locfdr"))'

## Install EpiDISH
RUN R -e 'BiocManager::install("EpiDISH")'

## Entry point
ENTRYPOINT [ "Rscript" ]