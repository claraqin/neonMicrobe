#' Set Base Directory for neonMicrobe Workspace
#'
#' Sets the specified directory as the workspace for the neonMicrobe pipeline. This base directory persists until the end of the R session. If data directories have not yet been generated within the base directory, use \code{\link{makeDataDirectories}} to do so.
#'
#' @param dir Directory to use as the base directory. Defaults to the working directory.
#'
#' @return No value is returned.
#' @export
#'
#' @seealso \code{\link{makeDataDirectories}}
#'
#' @examples
setBaseDirectory <- function(dir = getwd()) {
  if(dir.exists(dir)) {
    assign("NEONMICROBE_DIR_BASE", dir, envir = neonmicrobe_env)
  } else {
    message("Specified directory does not exist")
  }
}

#' Make Data Directories
#'
#' Creates the directories necessary for storing output data, rooted in the specified directory.
#'
#' @param check_location Default TRUE. Whether to prompt the user to confirm that the data directory structure should be created in the current working directory.
#'
#' @return No value is returned.
#' @export
#'
#' @seealso \code{\link{setBaseDirectory}}
#'
#' @examples
#' \dontrun{
#' makeDataDirectories()
#' }
makeDataDirectories <- function(check_location = TRUE) {
  if(identical(check_location, TRUE)) {
    resp <- readline(paste("Create new data directory structure in current working directory (", getwd(), ")? y/n: ", sep=" "))
    if(!(resp %in% c("y","Y"))) return(invisible(NULL))
  }

  # Get dynamic directory names
  seq_dir <- NEONMICROBE_DIR_SEQUENCE()
  seqmeta_dir <- NEONMICROBE_DIR_SEQMETA()
  soil_dir <- NEONMICROBE_DIR_SOIL()
  taxref_dir <- NEONMICROBE_DIR_TAXREF()
  outputs_dir <- NEONMICROBE_DIR_OUTPUTS()
  batches_dir <- NEONMICROBE_DIR_BATCHES()

  # if(is.null(seq_dir) | seq_dir == "") seq_dir <- file.path(getwd(), "data", "raw_sequence")
  # if(is.null(seqmeta_dir) | seqmeta_dir == "") seqmeta_dir <- file.path(getwd(), "data", "sequence_metadata")
  # if(is.null(soil_dir) | soil_dir == "") soil_dir <- file.path(getwd(), "data", "soil")
  # if(is.null(taxref_dir) | taxref_dir == "") taxref_dir <- file.path(getwd(), "data", "tax_ref")
  # if(is.null(outputs_dir) | outputs_dir == "") outputs_dir <- file.path(getwd(), "outputs")
  # if(is.null(outputs_dir) | outputs_dir == "") outputs_dir <- file.path(getwd(), "batch_outputs")

  message("Building output directories:")
  message("Using '", seq_dir, "' as sequence subdirectory.")
  message("Using '", seqmeta_dir, "' as sequence metadata subdirectory.")
  message("Using '", soil_dir, "' as soil data subdirectory.")
  message("Using '", taxref_dir, "' as taxonomic reference data subdirectory.")
  message("Using '", outputs_dir, "' as output data subdirectory.")
  message("Using '", batches_dir, "' as batches output data subdirectory.")

  createDirIfNotExist <- function(dir) {
    if(!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  }

  # If preset output directories for raw and QC'd sequence metadata do not exist, create them
  createDirIfNotExist(seqmeta_dir)
  raw_seqmeta_dir <- file.path(seqmeta_dir, "raw_metadata")
  qc_seqmeta_dir <- file.path(seqmeta_dir, "qc_metadata")
  createDirIfNotExist(raw_seqmeta_dir)
  createDirIfNotExist(qc_seqmeta_dir)

  # If preset directories for ITS and 16S data do not exist, create them
  createDirIfNotExist(seq_dir)
  seq_its_dir <- file.path(seq_dir, "ITS")
  seq_16s_dir <- file.path(seq_dir, "16S")
  createDirIfNotExist(file.path(seq_its_dir, "0_raw"))
  createDirIfNotExist(file.path(seq_16s_dir, "0_raw"))

  # Create intermediary directories for ITS and 16S data in the middle
  # of being processed
  createDirIfNotExist(outputs_dir)
  processing_its_dir <- file.path(outputs_dir, "mid_process", "ITS")
  processing_16s_dir <- file.path(outputs_dir, "mid_process", "16S")
  createDirIfNotExist(file.path(processing_its_dir, "1_filtN"))
  createDirIfNotExist(file.path(processing_its_dir, "2_trimmed"))
  createDirIfNotExist(file.path(processing_its_dir, "3_filtered"))
  createDirIfNotExist(file.path(processing_its_dir, "4_seqtabs"))
  createDirIfNotExist(file.path(processing_16s_dir, "1_trimmed"))
  createDirIfNotExist(file.path(processing_16s_dir, "2_filtered"))
  createDirIfNotExist(file.path(processing_16s_dir, "3_seqtabs"))

  # Also create directories for read-tracking tables
  read_tracking_its_dir <- file.path(outputs_dir, "track_reads", "ITS")
  read_tracking_16s_dir <- file.path(outputs_dir, "track_reads", "16S")
  createDirIfNotExist(read_tracking_its_dir)
  createDirIfNotExist(read_tracking_16s_dir)

  # If preset directory for soil data does not exist, create it
  createDirIfNotExist(soil_dir)

  # If preset directory for taxonomic reference data does not exist, create it
  createDirIfNotExist(taxref_dir)

  # Lastly, create directory for batch outputs. newBatch() will take care
  # of the rest
  createDirIfNotExist(batches_dir)
}
# Original parameters for makeDataDirectories:
# @param seq_dir,seqmeta_dir,soil_dir,taxref_dir,outputs_dir Absolute paths to directory where each type of data will be written (sequence data, sequence metadata, soil data, taxonomic reference data, processing outputs). Defaults to the following preset parameters, respectively: PRESET_OUTDIR_SEQUENCE, PRESET_OUTDIR_SEQMETA, PRESET_OUTDIR_SOIL, PRESET_OUTDIR_TAXREF, PRESET_OUTDIR_OUTPUTS.

###############################
# Dynamic directory names:

NEONMICROBE_DIR_BASE <- function() {
  if("NEONMICROBE_DIR_BASE" %in% ls(envir = neonmicrobe_env)) {
    get("NEONMICROBE_DIR_BASE", envir = neonmicrobe_env)
  } else {
    warning("Base directory for neonMicrobe (NEONMICROBE_DIR_BASE) has not been set. Use setBaseDirectory() to do so.")
    getwd()
  }
}
# for sequence data (fastq files)
NEONMICROBE_DIR_SEQUENCE <- function() {
  file.path(NEONMICROBE_DIR_BASE(), "data", "raw_sequence/")
}
# for sequence metadata
NEONMICROBE_DIR_SEQMETA <- function() {
  file.path(NEONMICROBE_DIR_BASE(), "data", "sequence_metadata/")
}
# for soil data
NEONMICROBE_DIR_SOIL <- function() {
  file.path(NEONMICROBE_DIR_BASE(), "data", "soil/")
}
# for taxonomy reference data
NEONMICROBE_DIR_TAXREF <- function() {
  file.path(NEONMICROBE_DIR_BASE(), "data", "tax_ref/")
}

# for outputs (sequence table, taxonomy table, phyloseq object). Contains mid_process and track_reads
NEONMICROBE_DIR_OUTPUTS <- function(dir = NULL) {
  checkArgsAgainstBatchParams("dir" = "DIR_OUTPUTS", verbose=FALSE)
  if(is.null(dir)) {
    file.path(NEONMICROBE_DIR_BASE(), "outputs")
  } else {
    dir
  }
}

# for batch-specific outputs (when processing batch is set using setBatch() or newBatch())
NEONMICROBE_DIR_BATCHES <- function() {
  file.path(NEONMICROBE_DIR_BASE(), "batch_outputs")
}

# NEONMICROBE_DIR_MIDPROCESS <- function(outputs_dir = NULL) { # for fastq files in the middle of processing into sequence tables
#   checkArgsAgainstBatchParams("outputs_dir" = "DIR_OUTPUTS", verbose=FALSE)
#   if(is.null(outputs_dir)) {
#     file.path(NEONMICROBE_DIR_OUTPUTS(), "mid_process")
#   } else {
#     file.path(outputs_dir, "mid_process")
#   }
# }
# NEONMICROBE_DIR_TRACKREADS <- function(outputs_dir = NULL) { # for read-tracking tables
#   checkArgsAgainstBatchParams("outputs_dir" = "DIR_OUTPUTS", verbose=FALSE)
#   if(is.null(outputs_dir)) {
#     file.path(NEONMICROBE_DIR_OUTPUTS(), "track_reads")
#   } else {
#     file.path(outputs_dir, "track_reads")
#   }
# }
###############################
