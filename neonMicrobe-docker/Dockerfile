# build off ubuntu 18.04 LTS
FROM rocker/r-base:4.0.4

# set env var with apps for  debian apt-get install
ENV BUILD_DEPS \
	dialog \
    apt-utils \
    tzdata
# update and install debian packages
RUN apt-get update && apt-get install -y $BUILD_DEPS

# install and configure miniconda3 for cutadapt
# need to update this for new absolute path in container
# install miniconda
ENV PATH="/home/rstudio/miniconda3/bin:${PATH}"
ARG PATH="/home/rstudio/miniconda3/bin:${PATH}"
RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /home/rstudio/.conda \
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
  && echo 'alias cutadapt = "conda run --prefix '/home/rstudio/miniconda3/envs/cutadaptenv' cutadapt"'


# use R script to install neonMicrobe
COPY install_packages.R install_packages.R
RUN R -e install_packages.R


