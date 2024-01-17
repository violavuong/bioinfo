# Build PRMeth
## Base image
FROM r-base:4.3.1

## Data/time of build
ARG buildDate

## Requirements needed to run PRMeth
ADD PRMeth.R /source/

######### METADATA #########
# Image information
LABEL about.image="R base image with required libraries for PRMeth"
LABEL image.authors="Viola Meixian Vuong"
LABEL image.authors.email="ViolaMeixian.Vuong@UGent.be"
## Update when new R and libraries version are available
LABEL image.maintainer="Maisa R Ferro dos Santos"
LABEL image.maintainer.email="MaisaRenata.FerrodosSantos@ugent.be"
LABEL base.image="r-base:4.3.1"
LABEL base_image.added_libraries="quadprog=1.5-8"
LABEL image.version="0.1"
## Date/time of build
LABEL buildDate=${buildDate}
# Software information
LABEL software="PRMeth for DNA methylation deconvolution in a partially-based reference framework"
LABEL software.version="No info available"
LABEL about.article="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04893-7"
LABEL about.documentation="https://github.com/hedingqin/PRMeth"
LABEL about.tags="DNA methylation, Deconvolution, Partial-reference"
# Docker command examples for running this image
LABEL docker.run01.help="docker run --rm prmeth:0.1 /source/PRMeth.R -h"
LABEL docker.run01.deconvolution="docker run --rm -v path-to-data-directory:/mnt -it prmeth:0.1 /source/PRMeth.R /mnt/path-to-ref /mnt/path-to-samples modality /mnt/path-to-outputs"
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
 procps \
 git

## Install PRMeth dependencies
RUN R -e 'install.packages("quadprog")'

## Clone PRMeth 
RUN git clone https://github.com/hedingqin/PRMeth

## Entry point
ENTRYPOINT [ "Rscript" ]