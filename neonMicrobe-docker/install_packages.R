#!/user/bin/env Rscript

install.packages(c("BiocManager", "neonUtilities", "vegan", "zetadiv", 
          "deblur", "decontam", "tidyverse", "R.utils"), repos = "https://cloud.r-project.org/", 
          dependencies = TRUE)

# set BiocManager Version
BiocManager::install(version = "3.11", ask = FALSE)

# install rhdf5 for neonUtilities
BiocManager::install("rhdf5")
# install dada2 with biocmanager version 3.11
BiocManager::install("dada2", version = "3.11")
BiocManager::install("DESeq2")
BiocManager::install("DECIPHER")
BiocManager::install("phyloseq")

# install geoNEON
devtools::install_github("NEONScience/NEON-geolocation/geoNEON", dependencies=TRUE)

#install neonMicrobe
devtools::install_github("claraqin/neonMicrobe")
