# Code to generate fastq files used in examples

library(plyr)
library(dplyr)
library(neonUtilities)
library(R.utils)
library(neonMicrobe)
library(ShortRead)
library(Biostrings)
library(dada2)

# Depends on the output of seqmeta_greatplains.R
seqmeta_greatplains <- read.csv("data-raw/seqmeta_greatplains.csv")

subset <- with(seqmeta_greatplains,
               union(grep("BMI_Plate23WellE7_16S_BDNB6", rawDataFileName),
                     grep("BMI_Plate37WellA12_16S_BJ8RK", rawDataFileName)))

downloadRawSequenceData(seqmeta_greatplains[subset,], outDir = "inst/extdata/", checkSize = FALSE)
