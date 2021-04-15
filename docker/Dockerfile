# build off ubuntu 18.04 LTS
FROM rocker/r-base:4.0.4

# set env var with apps for  debian apt-get install
ENV BUILD_DEPS \
	dialog \
    apt-utils \
    tzdata \
    vim \
    nano \
    curl \
    wget \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev 

# update and install debian packages
RUN apt-get update && apt-get install -y $BUILD_DEPS

# install and configure miniconda3 for cutadapt
# install miniconda
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

# update miniconda3
RUN conda init bash \
  && . ~/.bashrc \
  && conda update -n base -c defaults conda \
  && conda config --add channels defaults \
  && conda config --add channels bioconda \
  && conda config --add channels conda-forge \
  && conda install cutadapt \
  && echo 'alias cutadapt = "conda run --prefix '/miniconda3/envs/cutadaptenv' cutadapt"'


# use R exec to install neonMicrobe dependencies 
  # and other programs

RUN R -e 'install.packages(c("BiocManager", "tidyverse", "neonUtilities", "vegan", \
          "zetadiv", "testthat", "deblur", "decontam", "R.utils"), \ 
          repos = "https://cloud.r-project.org/", dependencies = TRUE)'

# set BiocManager Version
RUN R -e 'BiocManager::install(version = "3.11", ask = FALSE)'

# install rhdf5 for neonUtilities
RUN R -e 'BiocManager::install("rhdf5")'
# install dada2 with biocmanager version 3.11
RUN R -e 'BiocManager::install("dada2", version = "3.11")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("DECIPHER")'
RUN R -e 'BiocManager::install("phyloseq")'

#install devtools
RUN R -e 'install.packages("devtools", repos = "https://cloud.r-project.org/", dependencies = TRUE)'
# install geoNEON from GitHub with devtools
RUN R -e 'devtools::install_github("NEONScience/NEON-geolocation/geoNEON", dependencies=TRUE)'

#install neonMicrobe from GitHub with devtools
RUN R -e 'devtools::install_github("claraqin/neonMicrobe", dependencies=TRUE)'