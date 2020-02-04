# Handle raw sequence files downloaded by downloadRawSequenceData()

# This script automates the instructions in the "New server setup"
# page in the GitHub Wiki, EXCEPT that it only downloads ITS data
# (so far) and it does not perform the operations in the 
# "Pre-processing" section. 
# This script now splits files into ITS and 16S subdirectories.
# Furthermore, this script appends the sequencing run ID to all
# fastq filenames, so they can be processed by DADA2 in batches
# corresponding to sequencing runs.

source("./code/params.R")

BASE_DIR <- preset_outdir_sequence

meta <- downloadAllSequenceMetadata()
file_by_runid <- distinct(meta[,c("rawDataFileName","sequencerRunID")])

# get all the zip files
zipF <- list.files(path = BASE_DIR, pattern = "*.gz", full.names = TRUE)
zipF

# get the sequencer run ID associated with each zip file
runID_ind <- match(basename(zipF), file_by_runid[,"rawDataFileName"])
runIDs <- file_by_runid$sequencerRunID[runID_ind]
runIDs

# unzip files and append run ID to beginning of basenames
for(i in 1:length(zipF[2:6])) { # <-------------------------- TODO: CHANGE TO BE FOR ALL FILES
  # get run ID associated with the zip file
  runID <- runIDs[i]
  
  # unzip files
  untar(zipF[i], list=FALSE, exdir = BASE_DIR)
  
  # list all unzipped files
  # This regex matches all files with basenames that do NOT begin with "run" and which end with ".fastq"
  unzippedF_ITS <- grep("/(?!run)[^/]*_ITS[^/]*fastq$", list.files(path = BASE_DIR, recursive=TRUE, full.names=TRUE), perl=TRUE, value=TRUE)
  unzippedF_16S <- grep("/(?!run)[^/]*_16S[^/]*fastq$", list.files(path = BASE_DIR, recursive=TRUE, full.names=TRUE), perl=TRUE, value=TRUE)
  
  # rename all unzipped files by appending sequencer run ID
  unzippedF_ITS_rename <- file.path(BASE_DIR, "ITS", paste0("run", runID, "_", basename(unzippedF_ITS)))
  unzippedF_16S_rename <- file.path(BASE_DIR, "16S", paste0("run", runID, "_", basename(unzippedF_16S)))
  file.rename(unzippedF, unzippedF_rename)
}

list.files(path=file.path(BASE_DIR), recursive=TRUE)

# TODO: remove "hpc" directory

file.remove(list.files(path=BASE_DIR, pattern = "*.fastq$", recursive=TRUE, full.names=TRUE))
