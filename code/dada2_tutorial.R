library(dada2)

path <- "~/Downloads/R1"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = FALSE)) # If set to be FALSE, then working directory must contain the files
# fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = FALSE)) # none for now
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sub("_R1.fastq", "", fnFs)

plotQualityProfile(fnFs[1:2])
