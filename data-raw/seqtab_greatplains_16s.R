# Create collapsed NEON 16S sequence table for Great Plains sites

library(plyr)
library(dplyr)
library(neonUtilities)
library(R.utils)
library(neonMicrobe)
library(ShortRead)
library(Biostrings)
library(dada2)

meta_16s <- downloadSequenceMetadata(startYrMo = "2017-07", endYrMo = "2017-07",
                                     sites = c("KONZ", "CPER", "NOGP"),
                                     targetGene = "16S")

meta_16s_qc <- qcMetadata(meta_16s)

table(meta_16s_qc$siteID)
with(meta_16s_qc, table(siteID, sequencerRunID))

# download_success_16s <- downloadRawSequenceData(meta_16s_qc[which(meta_16s_qc$siteID=="KONZ"),])
download_success_16s <- downloadRawSequenceData(meta_16s_qc)

outdir_sequence <- PRESET_OUTDIR_SEQUENCE
fn0_16s <- file.path(outdir_sequence, meta_16s_qc$rawDataFileName)
fn_16s <- unique(fn0_16s[file.exists(fn0_16s)])

reorganized_files_16s <- organizeRawSequenceData(fn_16s, meta_16s_qc)


# Process sequences

if(is.null(PRESET_OUTDIR_SEQUENCE) | PRESET_OUTDIR_SEQUENCE == "") {
  PATH_RAW <- file.path(getwd(), "data", "raw_sequence", "16S", "0_raw")
} else {
  PATH_RAW <- file.path(PRESET_OUTDIR_SEQUENCE, "16S", "0_raw")
}

today <- Sys.Date()
if(is.null(PRESET_OUTDIR_OUTPUTS) | PRESET_OUTDIR_OUTPUTS == "") {
  PATH_OUTPUTS <- file.path(getwd(), "outputs", today)
} else {
  PATH_OUTPUTS <- file.path(PRESET_OUTDIR_OUTPUTS, today)
}

PATH_MIDPROCESS <- file.path(PATH_OUTPUTS, "mid_process", "16S")
PATH_TRACK <- file.path(PATH_OUTPUTS, "track_reads", "16S")

PATH_TRIMMED <- file.path(PATH_MIDPROCESS, "1_trimmed")
PATH_FILTERED <- file.path(PATH_MIDPROCESS, "2_filtered")
PATH_SEQTABS <- file.path(PATH_MIDPROCESS, "3_seqtabs")

dir.create(PATH_MIDPROCESS, recursive=TRUE)
dir.create(PATH_TRACK, recursive=TRUE)
dir.create(PATH_TRIMMED, recursive=TRUE)
dir.create(PATH_FILTERED, recursive=TRUE)
dir.create(PATH_SEQTABS, recursive=TRUE)

meta <- meta_16s_qc
meta_fn <- matchFastqToMetadata(list.files(PATH_RAW, full.names=TRUE), meta)

unique_runs <- unique(meta_fn$sequencerRunID)
unique_runs

seqtab_filenames <- c()
t1 <- Sys.time()
ti <- c()
for (i in 1:length(unique_runs)) {
  runID <- unique_runs[i]
  message(paste0("Began processing ", runID, " at ", Sys.time()))

  # Get fastq files from this sequencing run ID
  fn <- meta_fn$file[which(meta_fn$sequencerRunID==runID)]
  fn_base <- basename(fn)
  if(length(fn)==0) {
    message("No files found for sequencing run ", runID, ". Trying next sequencing run.")
    next
  }

  # Trim reads based on the primer lengths supplied in params.r
  trim_trackReads <- trimPrimers16S2(fn, PATH_TRIMMED, meta, "CCTACGGGNBGCASCAG", "GACTACNVGGGTATCTAATCC")

  # Filter reads based on the settings in params.r
  filter_trackReads <- qualityFilter16S2(file.path(PATH_TRIMMED, fn_base), PATH_FILTERED, meta, multithread = MULTITHREAD, maxEE = c(MAX_EE_FWD, MAX_EE_REV), truncLen=220)

  # Now create sequence table for run
  seqtab.list <- runDada16S2(file.path(PATH_FILTERED, fn_base), meta, remove_chimeras=TRUE, multithread=MULTITHREAD, verbose=VERBOSE)

  # Create output tracking file
  track <- Reduce(
    function(x, y, ...) transform(merge(x, y, by = 0, all = TRUE, ...), row.names=Row.names, Row.names = NULL),
    list(trim_trackReads,
         filter_trackReads[,2,drop=FALSE],
         seqtab.list$track)
  )
  names(track)[3] <- "filtered"
  track[is.na(track)] <- 0

  # Save tracking table (which tracks no. of reads remaining at each stage) and sequence table
  write.csv(track, file.path(PATH_TRACK, paste0("track_reads_", runID, ".csv")))
  saveRDS(seqtab.list$seqtab, file.path(PATH_SEQTABS, paste0("NEON_16S_seqtab_nochim_", runID, ".rds")))
  seqtab_filenames <- c(seqtab_filenames, file.path(PATH_SEQTABS, paste0("NEON_16S_seqtab_nochim_", runID, ".rds")))
  message(paste0("Finished processing reads in ", runID, " at ", Sys.time()))
  message(paste0("Sequencing run-specific sequence tables can be found in ", PATH_SEQTABS))

  # Clean up
  rm(trim_trackReads)
  rm(filter_trackReads)
  rm(seqtab.list)

  ti <- c(ti, Sys.time())
}

seqtab_joined <- mergeSequenceTables(tables = seqtab_filenames)

# Collapses NEON 16S ASV table

t0 <- Sys.time()
seqmeta_greatplains_16s <- collapseNoMismatch(seqtab_joined) # Took 6.93 hours on socs-stats.ucsc.edu
usethis::use_data(seqmeta_greatplains_16s, overwrite=TRUE)
t1 <- Sys.time()

t1 - t0
