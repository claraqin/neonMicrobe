# DADA2 tutorial
# Follows the ITS workflow https://benjjneb.github.io/dada2/ITS_workflow.html

library(dada2)

path <- "~/zhulab/NEON_DoB_analysis/data/Illumina/ITS"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = FALSE)) # If set to be FALSE, then working directory must contain the files
# fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = FALSE)) # none for now
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sub("_R1.fastq", "", fnFs)

plotQualityProfile(paste(path, fnFs[1:2], sep="/"))

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

