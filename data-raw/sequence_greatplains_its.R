# Code to generate fastq files used in examples

library(plyr)
library(dplyr)
library(neonUtilities)
library(R.utils)
library(neonMicrobe)
library(ShortRead)
library(Biostrings)
library(dada2)

# Depends on the output of seqmeta_greatplains_its.R
seqmeta_greatplains_its <- read.csv("data-raw/seqmeta_greatplains_its.csv")

subset <- with(seqmeta_greatplains_its,
               union(grep("BMI_Plate23WellG5_ITS_BMCC4", rawDataFileName),
                     grep("BMI_Plate37WellB2_ITS_BNM6G", rawDataFileName)))

downloadRawSequenceData(seqmeta_greatplains_its[subset,], outDir = "inst/extdata/", checkSize = FALSE)
