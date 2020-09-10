# DADA2 workflow for processing NEON ITS raw sequences
# Follows https://benjjneb.github.io/dada2/ITS_workflow.html

# This version of the ITS workflow process R1 only, and uses
# the new functions in ./R/utils.R

t0 <- Sys.time()

# Load parameters from params.R
source("./code/params.R")

# Load utility functions
source("./R/utils.R")

# Generate filepath names
if(is.null(PRESET_OUTDIR_SEQUENCE) | PRESET_OUTDIR_SEQUENCE == "") {
  PATH_ITS <- file.path(PRESET_OUTDIR, "raw_sequence", "ITS")
} else {
  PATH_ITS <- file.path(PRESET_OUTDIR, PRESET_OUTDIR_SEQUENCE, "ITS")
}
PATH_RAW <- file.path(PATH_ITS, "0_raw")
PATH_FILTN <- file.path(PATH_ITS, "1_filtN")
PATH_TRIMMED <- file.path(PATH_ITS, "2_trimmed")
PATH_FILTERED <- file.path(PATH_ITS, "3_filtered")
PATH_SEQTABS <- file.path(PATH_ITS, "4_seqtabs")
PATH_TRACK <- file.path(PATH_ITS, "track_reads")

# Load libraries
library(dada2)
library(ShortRead)
library(Biostrings)
library(tibble)
library(dplyr)

# Get all run IDs so you can group by them
unique_runs <- unique(unlist(
  regmatches(list.files(PATH_RAW), gregexpr("^run[A-Za-z0-9]*", list.files(PATH_RAW)))
))

# If SMALL_SUBSET == TRUE, run only the first runID
if (SMALL_SUBSET) {
  loop_length <- 1
} else {
  loop_length <- length(unique_runs)
}

all.seqtabs <- list()
t1 <- Sys.time()
ti <- c()
first <- TRUE
for (i in 1:loop_length) {
  runID <- unique_runs[i]
  message(paste0("Began processing ", runID, " at ", Sys.time()))

  # Forward fastq filenames have format: SAMPLENAME_R1_001.fastq
  fnFs <- sort(list.files(PATH_RAW, pattern=paste0(runID, ".*_R1.fastq"), full.names = TRUE))

  # If SMALL_SUBSET == TRUE,
  # keep only the first two forward-reverse pairs of sequence files
  if(SMALL_SUBSET){
    if(length(fnFs > 2)) fnFs <- fnFs[1:2]
  }

  fnFs_base <- basename(fnFs)

  # "Pre-filter" the sequences just to remove those with Ns, but perform no other filtering
  prefilter_trackReads <- qualityFilterITS(fnFs_base, PATH_RAW, PATH_FILTN, maxN = 0)

  # Trim primers from ITS sequences
  trimPrimersITS(fnFs_base, PATH_FILTN, PATH_TRIMMED, "CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC", discard_untrimmed = TRUE)

  # Filter and truncate
  filter_trackReads <- qualityFilterITS(fnFs_base, PATH_TRIMMED, PATH_FILTERED, MULTITHREAD, MAX_EE_FWD, TRUNC_Q, MIN_LEN)

  # Now create sequence table for run
  seqtab.list <- runDadaITS(fnFs_base, PATH_FILTERED, MULTITHREAD, VERBOSE)

  # Create output tracking file
  track <- cbind.data.frame(prefilter_trackReads,
                            filtered = filter_trackReads,
                            seqtab.list$track)
  names(track)[1:4] <- c("reads.in", "prefiltered.reads", "trimmed.reads", "filtered.reads")

  # Append sequence table to output list
  all.seqtabs[[runID]] <- seqtab.list$seqtab.nochim

  # Save tracking table (which tracks no. of reads remaining at each stage) and sequence table
  if(SMALL_SUBSET) {
    write.csv(track, file.path(PATH_TRACK, paste0("track_reads_",runID,"_SMALLSUBSET.csv")))
    saveRDS(seqtab.list$seqtab.nochim, file.path(PATH_SEQTABS, paste0("NEON_ITS_seqtab_nochim_", runID, "_SMALLSUBSET.rds")))
  } else {
    write.csv(track, file.path(PATH_TRACK, paste0("track_reads_",runID,".csv")))
    saveRDS(seqtab.list$seqtab.nochim, file.path(PATH_SEQTABS, paste0("NEON_ITS_seqtab_nochim_", runID, ".rds")))
  }
  message(paste0("Finished processing reads in ", runID, " at ", Sys.time()))
  message(paste0("Sequencing run-specific sequence tables can be found in ", PATH_SEQTABS))

  ti <- c(ti, Sys.time())
}
t2 <- Sys.time()

# Merge the sequence tables from all runs
seqtab_joined <- mergeSequenceTables(tables = all.seqtabs)
# Assign taxonomy using the UNITE reference database
tax <- assignTaxonomy(seqtab_joined, UNITE_REF_PATH, multithread = MULTITHREAD, tryRC = TRUE, verbose = VERBOSE)

# Save OTU table and taxonomic table as RDS files
# to hand off to dada2_to_phyloseq.R
saveRDS(seqtab_joined, file.path(PRESET_OUTDIR_DADA2, PRESET_FILENAME_JOINED_SEQTAB))
saveRDS(tax, file.path(PRESET_OUTDIR_DADA2, PRESET_FILENAME_TAXTAB))

