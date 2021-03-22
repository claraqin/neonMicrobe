#' Make Output Directories
#'
#' @param seq_dir,seqmeta_dir,soil_dir,taxref_dir,outputs_dir Absolute paths to directory where each type of data will be written (sequence data, sequence metadata, soil data, taxonomic reference data, processing outputs). Defaults to the following preset parameters, respectively: PRESET_OUTDIR_SEQUENCE, PRESET_OUTDIR_SEQMETA, PRESET_OUTDIR_SOIL, PRESET_OUTDIR_TAXREF, PRESET_OUTDIR_OUTPUTS.
#'
#' @return No value is returned.
#' @export
#'
#' @examples
#' \dontrun{
#' makeOutputDirectories()
#' }
makeOutputDirectories <- function(seq_dir=PRESET_OUTDIR_SEQUENCE,
                                  seqmeta_dir=PRESET_OUTDIR_SEQMETA, soil_dir=PRESET_OUTDIR_SOIL,
                                  taxref_dir=PRESET_OUTDIR_TAXREF, outputs_dir=PRESET_OUTDIR_OUTPUTS) {
  if(is.null(seq_dir) | seq_dir == "") seq_dir <- file.path(getwd(), "data", "raw_sequence")
  if(is.null(seqmeta_dir) | seqmeta_dir == "") seqmeta_dir <- file.path(getwd(), "data", "sequence_metadata")
  if(is.null(soil_dir) | soil_dir == "") soil_dir <- file.path(getwd(), "data", "soil")
  if(is.null(taxref_dir) | taxref_dir == "") taxref_dir <- file.path(getwd(), "data", "tax_ref")
  if(is.null(outputs_dir) | outputs_dir == "") outputs_dir <- file.path(getwd(), "outputs")

  message("Building output directories:")
  message("Using '", seq_dir, "' as sequence subdirectory.")
  message("Using '", seqmeta_dir, "' as sequence metadata subdirectory.")
  message("Using '", soil_dir, "' as soil data subdirectory.")
  message("Using '", taxref_dir, "' as taxonomic reference data subdirectory.")
  message("Using '", outputs_dir, "' as output data subdirectory.")

  createDirIfNotExist <- function(dir) {
    if(!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  }

  # If preset output directories do not exist, create them
  createDirIfNotExist(seq_dir)
  createDirIfNotExist(seqmeta_dir)
  createDirIfNotExist(soil_dir)
  createDirIfNotExist(taxref_dir)
  createDirIfNotExist(outputs_dir)

  # If preset output directories for raw and QC'd sequence metadata do not exist, create them
  raw_seqmeta_dir <- file.path(seqmeta_dir, "raw_metadata")
  qc_seqmeta_dir <- file.path(seqmeta_dir, "qc_metadata")
  createDirIfNotExist(raw_seqmeta_dir)
  createDirIfNotExist(qc_seqmeta_dir)

  # If preset directories for ITS and 16S data do not exist, create them
  seq_its_dir <- file.path(seq_dir, "ITS")
  seq_16s_dir <- file.path(seq_dir, "16S")
  createDirIfNotExist(file.path(seq_its_dir, "0_raw"))
  createDirIfNotExist(file.path(seq_16s_dir, "0_raw"))

  # Create intermediary directories for ITS and 16S data in the middle
  # of being processed
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
}
