#' QC Raw Sequence Metadata
#'
#' Performs basic QC checks and reports potential issues with downstream processing.
#'
#' @param metadata The output of downloadSequenceMetadata(). Must be provided as either the data.frame returned by downloadSequenceMetadata() or as a filepath to the csv file produced by downloadSequenceMetadata() when dir is specified.
#' @param outdir_sequence Default file.path(PRESET_OUTDIR, PRESET_OUTDIR_SEQUENCE). Directory where raw sequence files can be found before reorganizing.
#'
#' @return Report file of potential quality issues with metadata records
qcMetadata <- function(metadata, outdir = file.path(PRESET_OUTDIR, PRESET_OUTDIR_SEQMETA)) {

  metadata_load_err <- FALSE
  
  if(class(metadata) == "data.frame") {
    metadata <- metadata
  } else if(class(metadata) == "character") {
    if(file.exists(metadata)) {
      metadata <- read.csv(metadata)
    } else {
      metadata_load_err <- TRUE
    }
  } else {
    metadata_load_err <- TRUE
  }
  
  if(metadata_load_err) {
    stop("'metadata' must be the data.frame output from downloadSequenceMetadata() or a filepath to the local copy of the output from downloadSequenceMetadata()")
  }
  
  # QC check #1: completeness
  
  
  # Using metadata, look up the sequencing run ID and target gene for each file
  match_fn_to_meta <- match(basename(fn), metadata$rawDataFileName)
  runIDs <- metadata$sequencerRunID[match_fn_to_meta]
  targetGenes <- as.character(metadata$targetGene[match_fn_to_meta])
  targetGenes[grep("16S", targetGenes, ignore.case=TRUE)] <- "16S"
  targetGenes[grep("ITS", targetGenes, ignore.case=TRUE)] <- "ITS"
  
  if(all(is.na(runIDs))) stop("All sequencer run ID values are NA. Cannot proceed with reorganizing this download batch.
                              Try a different `fn` or `metadata`.")
  
  files_organized <- c()
  for(i in 1:length(fn)) {
    # get run ID associated with the zip file
    runID <- runIDs[i]
    if(is.na(runID)) {
      warning("No recorded sequencer run ID in metadata for file: ", fn[i], ". Skipping this file.")
      next
    }
    
    # get target gene associated with the zip file
    targetGene <- targetGenes[i]
    if(is.na(runID)) {
      warning("No recorded target gene in metadata for file: ", fn[i], ". Skipping this file.")
      next
    }
    if(!(targetGene %in% c("16S", "ITS"))) {
      warning("Recorded target gene in metadata for file does not clearly indicate '16S' or 'ITS': ", fn[i],
              ". Skipping this file.")
      next
    }
    
    # If the file is a tar file (unusual case)
    if(grepl("tar.gz", fn[i])) {
      tar_filenames <- untar(fn[i], list=TRUE)
      untarred_dir <- file.path(outdir_sequence, "untarred")
      untar_dir <- if(!dir.exists(untarred_dir)) dir.create(untarred_dir)
      untar(fn[i], list=FALSE, exdir = untarred_dir)
      file.remove(fn[i])
      untarred <- file.path(untarred_dir, tar_filenames)
      
      # rename untarred files by appending sequencer run ID and moving to target gene-specific subdirectory
      untarred_rename_to <- file.path(outdir_sequence, targetGene, "0_raw", paste0("run", runID, "_", basename(untarred)))
      if(length(untarred) > 0) {
        untarred_rename_to <- file.path(outdir_sequence, targetGene, "0_raw", paste0("run", runID, "_", basename(untarred)))
        renamed <- file.rename(untarred, untarred_rename_to)
        if(!all(renamed)) {
          warning(sum(!renamed), " untarred file(s) in sequencer run ", runID, "were not successfully reorganized")
        }
        files_organized <- c(files_organized, untarred_rename_to[renamed])
      }
      message("Extracted files from ", basename(fn[i]), " and reorganized its contents.")
      
      # If the file is not tarred (typical case)
    } else {
      
      # rename file by appending sequencer run ID and moving to target gene-specific subdirectory
      rename_to <- file.path(outdir_sequence, targetGene, "0_raw", paste0("run", runID, "_", basename(fn[i])))
      renamed <- file.rename(fn[i], rename_to)
      
      if(!renamed) {
        warning(fn[i], " was not successfully reorganized")
      } else {
        files_organized <- c(files_organized, rename_to)
      }
      message("Reorganized ", basename(fn[i]))
    }
  }
  
  # Gzip all non-gzipped fastq files
  for(i in 1:length(files_organized)) {
    if(!grepl(".gz$", files_organized[i])) {
      gzip(files_organized[i])
    }
  }
  
  if(length(files_organized) > 0) {
    message(length(files_organized), " out of ", length(fn), " files in sequencer run ", runID, " were successfully reorganized.")
    message("Reorganized files can be found in ", file.path(outdir_sequence, "ITS", "0_raw"), " and ", file.path(outdir_sequence, "16S", "0_raw"), ".")
    return(files_organized)
  } else {
    warning("No files were successfully reorganized")
    return(NULL)
  }
}
