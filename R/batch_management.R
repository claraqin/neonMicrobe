
#' Create New Processing Batch
#'
#'
#'
#' @param batch_id Character string. Unique ID to use for the new processing batch.
#' @param seq_meta_file Character string. File path to the sequence metadata to associate with this processing batch. Once set, this cannot be changed except by overwriting this batch.
#' @param batches_dir Default file.path(PRESET_OUTDIR_OUTPUS, "batches"). Directory where batch-specific directories are saved.
#' @param overwrite Default FALSE. If processing batch already exists in the specified directory, whether to overwrite it.
#' @param set_batch Default TRUE. Whether to switch to the new processing batch after creating it.
#'
#' @return No value returned
#' @export
#'
#' @examples
#' \dontrun{
#' newBatch("abc", "data/sequence_metadata/mmg_soilMetadata_ITS_2021-03-08134134.csv")
#' newBatch("xyz", "data/sequence_metadata/mmg_soilMetadata_ITS_2021-03-08134134.csv")
#' }
newBatch <- function(batch_id, seq_meta_file, batches_dir = file.path(PRESET_OUTDIR_OUTPUTS, "batches"), overwrite=FALSE, set_batch=TRUE) {
  if(!dir.exists(batches_dir)) dir.create(batches_dir, recursive=TRUE)

  metadata_load_err <- FALSE
  if(!file.exists(seq_meta_file)) {
    warning("Processing batch could not be created: No sequence metadata found at specified file location.")
  } else {
    this_batch_dir <- file.path(batches_dir, batch_id)
    if(dir.exists(this_batch_dir) & !identical(overwrite, TRUE)) {
      warning("Directory for processing batch ID ", batch_id, " already exists: ", this_batch_dir, ". To really overwrite, set overwrite=TRUE.")
    } else {
      # Overwrite if batch already exists and overwrite==TRUE
      if(dir.exists(this_batch_dir) & identical(overwrite, TRUE)) unlink(this_batch_dir, recursive=TRUE)
      # Create new batch
      dir.create(this_batch_dir)
      # Make ABOUT file in batch directory
      about_file <- file(file.path(this_batch_dir, "ABOUT"))
      writeLines(c(paste0("Batch created: ", Sys.time()),
                   paste0("Sequence metadata: ", normalizePath(seq_meta_file))),
                 about_file)
      close(about_file)

      # Make directory structure in new batch directory
      createDirIfNotExist <- function(dir) {
        if(!dir.exists(dir)) dir.create(dir, recursive=TRUE)
      }

      # First, create intermediary directories for ITS and 16S data in the middle
      # of being processed
      processing_its_dir <- file.path(this_batch_dir, "mid_process", "ITS")
      processing_16s_dir <- file.path(this_batch_dir, "mid_process", "16S")
      createDirIfNotExist(file.path(processing_its_dir, "1_filtN"))
      createDirIfNotExist(file.path(processing_its_dir, "2_trimmed"))
      createDirIfNotExist(file.path(processing_its_dir, "3_filtered"))
      createDirIfNotExist(file.path(processing_its_dir, "4_seqtabs"))
      createDirIfNotExist(file.path(processing_16s_dir, "1_trimmed"))
      createDirIfNotExist(file.path(processing_16s_dir, "2_filtered"))
      createDirIfNotExist(file.path(processing_16s_dir, "3_seqtabs"))

      # Also create directories for read-tracking tables
      read_tracking_its_dir <- file.path(this_batch_dir, "track_reads", "ITS")
      read_tracking_16s_dir <- file.path(this_batch_dir, "track_reads", "16S")
      createDirIfNotExist(read_tracking_its_dir)
      createDirIfNotExist(read_tracking_16s_dir)

      message("Created new batch directory at ", this_batch_dir)
      if(set_batch==TRUE) {
        setBatch(batch_id, batches_dir)
      }
    }
  }
}

#' Switch to Existing Processing Batch
#'
#' @param batch_id Character string. Unique ID of an existing processing batch.
#' @param batches_dir Default file.path(PRESET_OUTDIR_OUTPUS, "batches"). Directory where batch-specific directories are saved.
#'
#' @return
#' @export
#'
#' @examples
setBatch <- function(batch_id = NULL, batches_dir = file.path(PRESET_OUTDIR_OUTPUTS, "batches")) {
  # If no batch ID is provided, clear the WORKING_BATCH_ID variable
  if(is.null(batch_id) | missing(batch_id)) {
    if ("WORKING_BATCH_ID" %in% ls(envir = .GlobalEnv)) {
      prev_batch_id <- get("WORKING_BATCH_ID", envir = .GlobalEnv)
      rm("WORKING_BATCH_ID", envir = .GlobalEnv)
      rm("WORKING_BATCH_DIR", envir = .GlobalEnv)
      message("Stepped outside of processing batch structure. No working batch ID is currently set. To undo this, run setBatch('", prev_batch_id, "')")
    }

  # Otherwise, set WORKING_BATCH_ID to batch_id
  } else {
    if(!identical("character", class(batch_id))) {
      stop("Value provided was not a character string. Unable to set batch.")
    }
    if(!dir.exists(batches_dir)) dir.create(batches_dir, recursive=TRUE)
    this_batch_dir <- file.path(batches_dir, batch_id)
    if(dir.exists(this_batch_dir)) {
      message("Switched to pre-existing batch at", this_batch_dir)
    } else {
      message("Batch does not exist. Create a new batch using newBatch().")
    }
    WORKING_BATCH_ID <<- batch_id
    WORKING_BATCH_DIR <<- normalizePath(this_batch_dir)
    message("Now working with processing batch ID ", batch_id)
  }
}


#' Get Current Batch ID
#'
#' Get the unique ID of the current processing batch.
#'
#' @param verbose Default FALSE. If TRUE, returns information about the current processing batch.
#'
#' @return Character. Unique ID of the current processing batch.
#' @export
#'
#' @examples
#' getBatch() # "xyz"
#' getBatch(verbose = TRUE) # "xyz", and table containing info about batch "xyz"
getBatch <- function(verbose=FALSE) {
  if ("WORKING_BATCH_ID" %in% ls(envir = .GlobalEnv)) {
    batch_id <- get("WORKING_BATCH_ID", envir = .GlobalEnv)
    if(verbose==TRUE) {
      this_batch_dir <- get("WORKING_BATCH_DIR", envir = .GlobalEnv)
      about <- readLines(con = file.path(this_batch_dir, "ABOUT"), 2)
      about_parsed <- sub(".*: ", "", about)
      about_parsed[2] <- basename(about_parsed[2])
      names(about_parsed) <- c("created", "sequence_metadata")
      print(about_parsed)
    }
    return(batch_id)
  } else {
    message("Processing batch is not currently set. Use setBatch() to do this.")
    return(invisible(NULL))
  }
}

#' List Processing Batches
#'
#' @param batches_dir
#'
#' @return
#' @export
#'
#' @examples
listBatches <- function(batches_dir = file.path(PRESET_OUTDIR_OUTPUTS, "batches")) {
  if(dir.exists(batches_dir)) {
    message("List of processing batches in ", batches_dir, ":")
    dirs <- list.dirs(batches_dir, full.names=TRUE, recursive=FALSE)
    about <- t(sapply(dirs, function(x) readLines(con = file.path(x, "ABOUT"), 2)))
    about_parsed <- sub(".*: ", "", about)
    about_parsed[,2] <- basename(about_parsed[,2])
    return(data.frame(
      ID = basename(dirs),
      created = about_parsed[,1],
      sequence_metadata = about_parsed[,2],
      row.names = NULL
    ))
  } else {
    message("Specified directory does not exist. Create your first processing batch using newBatch().")
  }
}



#' Check Function Arguments against Batch-Specific Parameters
#'
#' If a processing batch is currently set, inserting this function into another function from
#' neonMicrobe will check the parent function's arguments against the parameters associated
#' with the current processing batch. If override=TRUE, the current batch's parameters would
#' take precedence. For example, rather than writing to the standard outputs directory, the
#' parent function would write to the current batch's outputs directory. Rather than using
#' quality filtering parameters defined on-the-fly, the parent function would use the quality
#' filtering parameters associated with the current batch.
#'
#' @return
#'
#' @examples
checkArgsAgainstBatchParams <- function(override=FALSE, ...) {
  # This might be useful https://stackoverflow.com/questions/42873592/assign-variable-in-parent-environment-of-a-function
  # Dots (...) will refer to any arguments that have values defined in parameters, set equal to whatever
  # value they have in the parent function, e.g. "MAX_EE_FWD = maxEE"
}
