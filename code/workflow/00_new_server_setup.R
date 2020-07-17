# Handle raw sequence files downloaded by downloadRawSequenceData()

# This script automates the instructions in the "New server setup"
# page in the GitHub Wiki, EXCEPT that it does not perform the
# operations in the "Pre-processing" section.
# This script now splits files into ITS and 16S subdirectories.
# Furthermore, this script appends the sequencing run ID to all
# fastq filenames, so they can be processed by DADA2 in batches
# corresponding to sequencing runs.

source("./code/params.R")
source("./code/utils.R")

BASE_DIR <- PRESET_OUTDIR_SEQUENCE

# If preset output directories do not exist, create them
if(!dir.exists(PRESET_OUTDIR_SEQUENCE)) dir.create(PRESET_OUTDIR_SEQUENCE, recursive=TRUE)
if(!dir.exists(PRESET_OUTDIR_SEQMETA)) dir.create(PRESET_OUTDIR_SEQMETA, recursive=TRUE)
if(!dir.exists(PRESET_OUTDIR_SOIL)) dir.create(PRESET_OUTDIR_SOIL, recursive=TRUE)

# If preset output directories for ITS and 16S data do not exist, create them
if(!dir.exists(file.path(BASE_DIR, "ITS"))) dir.create(file.path(BASE_DIR, "ITS"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "16S"))) dir.create(file.path(BASE_DIR, "16S"), recursive=TRUE)

# Create intermediary directories for ITS and 16S data in the middle
# of being processed
if(!dir.exists(file.path(BASE_DIR, "ITS", "0_unzipped"))) dir.create(file.path(BASE_DIR, "ITS", "0_unzipped"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "ITS", "1_filtN"))) dir.create(file.path(BASE_DIR, "ITS", "1_filtN"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "ITS", "2_cutadapt"))) dir.create(file.path(BASE_DIR, "ITS", "2_cutadapt"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "ITS", "3_filtered"))) dir.create(file.path(BASE_DIR, "ITS", "3_filtered"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "ITS", "4_seqtabs"))) dir.create(file.path(BASE_DIR, "ITS", "4_seqtabs"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "ITS", "track_reads"))) dir.create(file.path(BASE_DIR, "ITS", "track_reads"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "16S", "0_unzipped"))) dir.create(file.path(BASE_DIR, "16S", "0_unzipped"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "16S", "1_filtN"))) dir.create(file.path(BASE_DIR, "16S", "1_filtN"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "16S", "2_cutadapt"))) dir.create(file.path(BASE_DIR, "16S", "2_cutadapt"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "16S", "3_filtered"))) dir.create(file.path(BASE_DIR, "16S", "3_filtered"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "16S", "4_seqtabs"))) dir.create(file.path(BASE_DIR, "16S", "4_seqtabs"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "16S", "track_reads"))) dir.create(file.path(BASE_DIR, "16S", "track_reads"), recursive=TRUE)

# Download sequence data and metadata
meta <- downloadSequenceMetadata()
file_by_runid <- distinct(meta[,c("rawDataFileName","sequencerRunID")])
downloadRawSequenceData()

# get all the zip files
zipF <- list.files(path = BASE_DIR, pattern = "*.gz", full.names = TRUE)
zipF

# get the sequencer run ID associated with each zip file
runID_ind <- match(basename(zipF), file_by_runid[,"rawDataFileName"])
runIDs <- file_by_runid$sequencerRunID[runID_ind]

# unzip files and append run ID to beginning of basenames
for(i in 1:length(zipF)) { # <------------------------ TODO: Replace with foreach or map or other parallel process
  # get run ID associated with the zip file
  runID <- runIDs[i]

  # unzip files
  untar(zipF[i], list=FALSE, exdir = BASE_DIR)

  # list all unzipped files
  # This regex matches all files with basenames that do NOT begin with "run" and which end with ".fastq"
  unzippedF_ITS <- grep("/(?!run)[^/]*_ITS[^/]*fastq$", list.files(path = BASE_DIR, recursive=TRUE, full.names=TRUE), perl=TRUE, value=TRUE)
  unzippedF_16S <- grep("/(?!run)[^/]*_16S[^/]*fastq$", list.files(path = BASE_DIR, recursive=TRUE, full.names=TRUE), perl=TRUE, value=TRUE)

  # rename all unzipped files by appending sequencer run ID
  unzippedF_ITS_rename <- file.path(BASE_DIR, "ITS", "0_unzipped", paste0("run", runID, "_", basename(unzippedF_ITS)))
  unzippedF_16S_rename <- file.path(BASE_DIR, "16S", "0_unzipped", paste0("run", runID, "_", basename(unzippedF_16S)))
  if(length(unzippedF_ITS) > 0) file.rename(unzippedF_ITS, unzippedF_ITS_rename)
  if(length(unzippedF_16S) > 0) file.rename(unzippedF_16S, unzippedF_16S_rename)
}

# Remove "hpc" directory, where files were moved out of
unlink(file.path(BASE_DIR, "hpc"), recursive=TRUE)

# Need to clear all files?
# list.files(path=file.path(BASE_DIR), recursive=TRUE)
# file.remove(list.files(path=BASE_DIR, pattern = "*.fastq$", recursive=TRUE, full.names=TRUE))
