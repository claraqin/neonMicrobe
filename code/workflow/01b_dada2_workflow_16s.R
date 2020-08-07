# DADA2 workflow for processing NEON 16S raw sequences
# Follows https://benjjneb.github.io/dada2/tutorial.html

t0 <- Sys.time()

# Load parameters from params.R
source("./code/params.R")

# Load utility functions
# source("./code/utils.R")
# source("./code/utils_16S.r")
source("./R/utils.R")

# Generate filepath names
PATH_16S <- file.path(PRESET_OUTDIR_SEQUENCE, "16S")
PATH_UNZIPPED <- file.path(PATH_16S, "0_unzipped")
PATH_CUT <- file.path(PATH_16S, "1_trimmed")
PATH_FILTERED <- file.path(PATH_16S, "2_filtered")
PATH_SEQTABS <- file.path(PATH_16S, "3_seqtabs")
PATH_TRACK <- file.path(PATH_16S, "track_reads")

# Load libraries
library(dada2)
library(ShortRead)
library(Biostrings)

# Get all run IDs so you can group by them
unique_runs <- unique(unlist(
  regmatches(list.files(PATH_UNZIPPED), gregexpr("^run[A-Za-z0-9]*", list.files(PATH_UNZIPPED)))
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

  # Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
  fnFs <- sort(list.files(PATH_UNZIPPED, pattern=paste0(runID, ".*_R1.fastq"), full.names = TRUE))
  fnRs <- sort(list.files(PATH_UNZIPPED, pattern=paste0(runID, ".*_R2.fastq"), full.names = TRUE))


  # If SMALL_SUBSET == TRUE,
  # keep only the first two forward-reverse pairs of sequence files
  if(SMALL_SUBSET){
    if(length(fnFs > 2)) fnFs <- fnFs[1:2]
    if(length(fnRs > 2)) fnRs <- fnRs[1:2]
  }

  # Remove any files that only have forward or reverse reads
  matched_fn <- remove_unmatched_files(fnFs, fnRs)
  fnFs <- matched_fn[[1]]
  fnRs <- matched_fn[[2]]

  # Trim reads based on the primer lengths supplied in params.r
  trim_trackReads <- trimPrimers16S(fnFs, fnRs, PATH_CUT, PRIMER_16S_FWD, PRIMER_16S_REV, MULTITHREAD)

  # TODO: rather than taking input directory, may be necessary to take input files because
  #       we might not be interested in analyzing all files in the directory (e.g. those
  #       preexisting from previous processing batches).

  # Filter reads based on the settings in params.r
  filter_trackReads <- qualityFilter16S(PATH_CUT, PATH_FILTERED, MULTITHREAD, MAX_EE_FWD, MAX_EE_REV, TRUNC.LENGTHS = c(265, 210))

  # Now create sequence table for run
  seqtab.list <- runDada16S(PATH_FILTERED, MULTITHREAD, VERBOSE)

  # Create output tracking file
  track <- cbind.data.frame(trim_trackReads,
                            filtered = filter_trackReads[,2],
                            seqtab.list$track)

  # Append sequence table to output list
  all.seqtabs[[runID]] <- seqtab.list$seqtab.nochim

  # Save tracking table (which tracks no. of reads remaining at each stage) and sequence table
  if(SMALL_SUBSET) {
    write.csv(track, file.path(PATH_TRACK, paste0("track_reads_",runID,"_SMALLSUBSET.csv")))
    saveRDS(seqtab.list$seqtab.nochim, file.path(PATH_SEQTABS, paste0("NEON_16S_seqtab_nochim_", runID, "_SMALLSUBSET.rds")))
  } else {
    write.csv(track, file.path(PATH_TRACK, paste0("track_reads_",runID,".csv")))
    saveRDS(seqtab.list$seqtab.nochim, file.path(PATH_SEQTABS, paste0("NEON_16S_seqtab_nochim_", runID, ".rds")))
  }
  message(paste0("Finished tracking reads through pipeline in ", runID, " at ", Sys.time()))
  message(paste0("Finished saving sequence table of ", runID, " at ", Sys.time()))
  message(paste0("Sequencing run-specific sequence tables can be found in ", PATH_SEQTABS))
}

# Merge the sequence tables from all runs
seqtab_joined <- mergeSequenceTables(tables = all.seqtabs)
# Assign taxonomy using the Silva reference database
tax <- assignTaxonomy(seqtab_joined, SILVA_REF_PATH, multithread = MULTITHREAD, verbose = VERBOSE)

# Save OTU table and taxonomic table as RDS files
# to hand off to dada2_to_phyloseq.R
saveRDS(seqtab_joined, "./data/NEON_16S_seqtab_nochim.Rds")
saveRDS(tax, "./data/NEON_16S_tax.Rds")
