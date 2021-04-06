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
#' setBaseDirectory() # sets current working directory as base directory
#' setBaseDirectory(dirname(getwd())) # sets parent directory as base directory
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

  message("Building output directories:")
  message("Using '", seq_dir, "' as sequence subdirectory.")
  message("Using '", seqmeta_dir, "' as sequence metadata subdirectory.")
  message("Using '", soil_dir, "' as soil data subdirectory.")
  message("Using '", taxref_dir, "' as taxonomic reference data subdirectory.")
  message("Using '", outputs_dir, "' as output data subdirectory.")
  message("Using '", batches_dir, "' as batches output data subdirectory.")

  createDirIfNotExist <- function(dir) {
    if(length(dir) == 1) {
      if(!dir.exists(dir)) dir.create(dir, recursive=TRUE)
    } else {
      for(d in dir) {
        if(!dir.exists(d)) dir.create(d, recursive=TRUE)
      }
    }
  }

  # If preset output directories for raw and QC'd sequence metadata do not exist, create them
  createDirIfNotExist(file.path(seqmeta_dir, c("raw_metadata",
                                               "qc_metadata")))

  # If preset directories for ITS and 16S data do not exist, create them
  createDirIfNotExist(file.path(seq_dir, c("ITS",
                                           "16S")))

  # Create intermediary directories for ITS and 16S data in the middle
  # of being processed
  createDirIfNotExist(file.path(outputs_dir, "mid_process", "ITS", c("1_filtN",
                                                                     "2_trimmed",
                                                                     "3_filtered",
                                                                     "4_seqtabs")))
  createDirIfNotExist(file.path(outputs_dir, "mid_process", "16S", c("1_trimmed",
                                                                     "2_filtered",
                                                                     "3_seqtabs")))

  # Also create directories for read-tracking tables
  createDirIfNotExist(file.path(outputs_dir, "track_reads", c("ITS",
                                                              "16S")))

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

#' Dynamic Directory Name for Base Directory
#'
#' @return Directory path (character).
#' @export
NEONMICROBE_DIR_BASE <- function() {
  if("NEONMICROBE_DIR_BASE" %in% ls(envir = neonmicrobe_env)) {
    get("NEONMICROBE_DIR_BASE", envir = neonmicrobe_env)
  } else {
    warning("Base directory for neonMicrobe (NEONMICROBE_DIR_BASE) has not been set. Use setBaseDirectory() to do so.")
    getwd()
  }
}

#' Dynamic Directory Name for Raw Sequence Data
#'
#' For raw sequence data (fastq files) from NEON soil microbe marker gene sequences
#' (DP1.10108.001).
#'
#' @return Directory path (character).
#' @export
NEONMICROBE_DIR_SEQUENCE <- function() {
  file.path(NEONMICROBE_DIR_BASE(), "data", "raw_sequence/")
}

#' Dynamic Directory Name for Sequence Metadata
#'
#' For metadata associated with NEON soil microbe marker gene sequences (DP1.10108.001).
#'
#' @return Directory path (character).
#' @export
NEONMICROBE_DIR_SEQMETA <- function() {
  file.path(NEONMICROBE_DIR_BASE(), "data", "sequence_metadata/")
}

#' Dynamic Directory Name for Soil Data
#'
#' For NEON soil data DP1.10086.001:
#' "Soil physical and chemical properties, periodic",
#' tables sls_soilCoreCollection, sls_soilMoisture, sls_soilpH, and sls_soilChemistry.
#'
#' @return Directory path (character).
#' @export
NEONMICROBE_DIR_SOIL <- function() {
  file.path(NEONMICROBE_DIR_BASE(), "data", "soil/")
}


#' Dynamic Directory Name for Taxonomic Reference Data
#'
#' For taxonomic reference data, e.g. UNITE and SILVA databases.
#' Not strictly necessary for functioning, so long as user specifies
#' reference datasets in \code{\link[dada2]{assignTaxonomy}}.
#'
#' @return Directory path (character).
#' @export
NEONMICROBE_DIR_TAXREF <- function() {
  file.path(NEONMICROBE_DIR_BASE(), "data", "tax_ref/")
}

#' Dynamic Directory Name for Outputs
#'
#' For outputs (sequence table, taxonomy table, phyloseq object).
#' Contains mid_process and track_reads. Subject to modification
#' by batch processing parameters.
#'
#' @param dir If not NULL, directory path to force this to return.
#'
#' @return Directory path (character).
#' @export
NEONMICROBE_DIR_OUTPUTS <- function(dir = NULL) {
  checkArgsAgainstBatchParams("dir" = "DIR_OUTPUTS", verbose=FALSE)
  if(is.null(dir)) {
    file.path(NEONMICROBE_DIR_BASE(), "outputs")
  } else {
    dir
  }
}

#' Dynamic Directory Name for Mid-Processing Fastq Files
#'
#' For sequence files in the middle of being processed. Subject
#' to modification by batch processing parameters via their
#' modifications to \code{link{NEONMICROBE_DIR_OUTPUTS}}, of which
#' this is a subdirectory. Contains subdirs "16S", "ITS".
#'
#' @return Directory path (character).
#' @export
NEONMICROBE_DIR_MIDPROCESS <- function() {
  file.path(NEONMICROBE_DIR_OUTPUTS(), "mid_process")
}

#' Dynamic Directory Name for Processing Batch Outputs
#'
#' For batch-specific outputs (when processing batch is set using \code{\link{setBatch}} or \code{\link{newBatch}})
#'
#' @return Directory path (character).
#' @export
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

#' Dynamic Directory Name for Read-Tracking Tables
#'
#' For read-tracking tables. Subject to modification by batch
#' processing parameters via their modifications to \code{link{NEONMICROBE_DIR_OUTPUTS}},
#' of which this is a subdirectory. Contains subdirs "16S", "ITS".
#'
#' @return Directory path (character).
#' @export
NEONMICROBE_DIR_TRACKREADS <- function() { # for read-tracking tables
  file.path(NEONMICROBE_DIR_OUTPUTS(), "track_reads")
}
###############################
