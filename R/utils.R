# Utilities for use in other R scripts in NEON/DoB microbial community analysis

# This version of utils.R began as a copy of ./code/utils.R, and is intended for
# writing documentation of new functions using Roxygen. It does not include all
# of the functions in the original ./code/utils.R, like get.sample.name().
# Those functions are either copied from or closely adapted from the DADA tutorials,
# so they should be included in the vignettes, but not claimed as a new function.


#' Make Output Directories
#'
#' @param base_dir Default PRESET_OUTDIR in params.R. Absolute path to directory within which all output subdirectories will be nested.
#' @param seq_dirname Default PRESET_OUTDIR_SEQUENCE in params.R. Path to append to base_dir where sequence data will be output.
#' @param seqmeta_dirname Default PRESET_OUTDIR_SEQMETA in params.R. Path to append to base_dir where sequence metadata will be output.
#' @param soil_dirname Default PRESET_OUTDIR_SOIL in params.R. Path to append to base_dir where soil data will be output.
#'
#' @return No value is returned.
#'
#' @examples
#' \dontrun{
#' makeOutputDirectories(base_dir="/data/NEON")
#' }
makeOutputDirectories <- function(base_dir = PRESET_OUTDIR, seq_dirname=PRESET_OUTDIR_SEQUENCE,
                                  seqmeta_dirname=PRESET_OUTDIR_SEQMETA, soil_dirname=PRESET_OUTDIR_SOIL) {
  if(is.null(seq_dirname) | seq_dirname == "") seq_dirname <- "raw_sequence"
  if(is.null(seqmeta_dirname) | seqmeta_dirname == "") seqmeta_dirname <- "sequence_metadata"
  if(is.null(soil_dirname) | soil_dirname == "") soil_dirname <- "soil"

  seq_dir <- file.path(base_dir, seq_dirname)
  seqmeta_dir <- file.path(base_dir, seqmeta_dirname)
  soil_dir <- file.path(base_dir, soil_dirname)

  message("Building output directories from base directory '", base_dir, "'.")
  message("Using '", seq_dir, "' as sequence subdirectory.")
  message("Using '", seqmeta_dir, "' as sequence metadata subdirectory.")
  message("Using '", soil_dir, "' as soil data subdirectory.")

  seq_its_dir <- file.path(seq_dir, "ITS")
  seq_16s_dir <- file.path(seq_dir, "16S")

  # If preset output directories do not exist, create them
  if(!dir.exists(seq_dir)) dir.create(seq_dir, recursive=TRUE)
  if(!dir.exists(seqmeta_dir)) dir.create(seqmeta_dir, recursive=TRUE)
  if(!dir.exists(soil_dir)) dir.create(soil_dir, recursive=TRUE)

  # If preset output directories for ITS and 16S data do not exist, create them
  if(!dir.exists(seq_its_dir)) dir.create(seq_its_dir)
  if(!dir.exists(seq_16s_dir)) dir.create(seq_16s_dir)

  # Create intermediary directories for ITS and 16S data in the middle
  # of being processed
  if(!dir.exists(file.path(seq_its_dir, "0_raw"))) dir.create(file.path(seq_its_dir, "0_raw"))
  if(!dir.exists(file.path(seq_its_dir, "1_filtN"))) dir.create(file.path(seq_its_dir, "1_filtN"))
  if(!dir.exists(file.path(seq_its_dir, "2_trimmed"))) dir.create(file.path(seq_its_dir, "2_trimmed"))
  if(!dir.exists(file.path(seq_its_dir, "3_filtered"))) dir.create(file.path(seq_its_dir, "3_filtered"))
  if(!dir.exists(file.path(seq_its_dir, "4_seqtabs"))) dir.create(file.path(seq_its_dir, "4_seqtabs"))
  if(!dir.exists(file.path(seq_its_dir, "track_reads"))) dir.create(file.path(seq_its_dir, "track_reads"))
  if(!dir.exists(file.path(seq_16s_dir, "0_raw"))) dir.create(file.path(seq_16s_dir, "0_raw"))
  if(!dir.exists(file.path(seq_16s_dir, "1_trimmed"))) dir.create(file.path(seq_16s_dir, "1_trimmed"))
  if(!dir.exists(file.path(seq_16s_dir, "2_filtered"))) dir.create(file.path(seq_16s_dir, "2_filtered"))
  if(!dir.exists(file.path(seq_16s_dir, "3_seqtabs"))) dir.create(file.path(seq_16s_dir, "3_seqtabs"))
  if(!dir.exists(file.path(seq_16s_dir, "track_reads"))) dir.create(file.path(seq_16s_dir, "track_reads"))
}


#' Download NEON Marker Gene Sequencing Raw Data
#'
#' Downloads NEON raw sequence data files to the specified filepath, by referencing
#' URLs in the metadata output from \code{\link{downloadSequenceMetadata}}.
#'
#' @param metadata The output of downloadSequenceMetadata(). Must be provided as either the data.frame returned by downloadSequenceMetadata() or as a filepath to the csv file produced by downloadSequenceMetadata() when outdir is provided.
#' @param outdir Location where output files are saved. Defaults to file.path(PRESET_OUTDIR, PRESET_OUTDIR_SEQUENCE) in params.R.
#' @param ignore_tar_files If TRUE (default), does not download tar files. Each tar file is a batch containing an entire sequence run of fastq files. The tar file structure will soon be deprecated.
#' @param verbose If TRUE, prints status messages and progress bars associated with file downloads.
#'
#' @return Returns (invisibly) an integer code: 0 indicates success of downloads and a non-zero integer indicates failure. See the help page for \code{\link[utils]{download.file}} for more details.
downloadRawSequenceData <- function(metadata, outdir = file.path(PRESET_OUTDIR, PRESET_OUTDIR_SEQUENCE), ignore_tar_files=TRUE, verbose=FALSE) {

  library(utils)
  library(dplyr)

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

  if(ignore_tar_files) {
    tar_ind <- grep('\\.tar\\.gz', metadata$rawDataFileName)
    if(length(tar_ind) > 0) {
      metadata <- metadata[-tar_ind, ]
      message(length(tar_ind), " row(s) in metadata associated with batch-level sequence data
              were ignored prior to downloading raw sequence data.")
    }
  }

  u.urls <- unique(metadata$rawDataFilePath)
  fileNms <- gsub('^.*\\/', "", u.urls)
  print(paste("There are", length(u.urls), "unique (zipped) raw sequence files to download.") )

  for(i in 1:length(u.urls)) {
    download_success <- download.file(
      url = as.character(u.urls[i]),
      destfile = ifelse(dir.exists(outdir), paste(outdir, fileNms[i], sep="/"), paste(getwd(), fileNms[i], sep="/" )),
      quiet = !verbose
    )
    if(dir.exists(outdir)) {
      message("Finished downloading ", paste(outdir, fileNms[i], sep="/"))
    } else {
      message("Finished downloading ", paste(getwd(), fileNms[i], sep="/" ))
    }
  }

  return(download_success)
}


#' Organize Raw Sequence Data
#'
#' Moves raw sequence data into the correct subdirectory for the processing pipeline,
#' renames files to include sequencer run ID, and untars sequence data if necessary.
#'
#' @param fn Character vector of full names (including path) of raw sequence files. Can include tarballs.
#' @param metadata The output of downloadSequenceMetadata(). Must be provided as either the data.frame returned by downloadSequenceMetadata() or as a filepath to the csv file produced by downloadSequenceMetadata() when outdir is provided.
#' @param outdir_sequence Default file.path(PRESET_OUTDIR, PRESET_OUTDIR_SEQUENCE). Directory where raw sequence files are output.
#'
#' @return Character vector of the files (including files within tarballs) that were successfully reorganized. If no files were successfully reorganized, returns no value.
organizeRawSequenceData <- function(fn, metadata, outdir_sequence = file.path(PRESET_OUTDIR, PRESET_OUTDIR_SEQUENCE)) {
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

  runIDs <- metadata$sequencerRunID[match(basename(fn), metadata$rawDataFileName)]

  if(all(is.na(runIDs))) stop("All values of metadata$sequencerRunID are NA. Cannot append runID to beginning of filenames.")

  files_organized <- c()
  for(i in 1:length(fn)) {
    # get run ID associated with the zip file
    runID <- runIDs[i]
    if(is.na(runID)) next

    # If the file is a tar file (unusual case)
    if(grepl("tar.gz", fn[i])) {
      tar_filenames <- untar(fn[i], list=TRUE)
      untar(fn[i], list=FALSE, exdir = outdir_sequence)
      untarred <- file.path(outdir_sequence, tar_filenames)
      untarred_ITS <- grep("_ITS[^/]*fastq", untarred, value=TRUE)
      untarred_16S <- grep("_16S[^/]*fastq", untarred, value=TRUE)
      if(length(untarred_ITS) == 0 & length(untarred_16S) == 0) {
        warning("Could not distinguish between ITS and 16S files contained within ", fn[i], ". We use regex patterns '_ITS[^/]*fastq' and '_16S[^/]*fastq'
              on the untarred filenames to make this distinction.")
        next
      }

      # rename untarred files by appending sequencer run ID and moving to target gene-specific subdirectory
      untarred_ITS_rename_to <- file.path(outdir_sequence, "ITS", "0_raw", paste0("run", runID, "_", basename(untarred_ITS)))
      untarred_16S_rename_to <- file.path(outdir_sequence, "16S", "0_raw", paste0("run", runID, "_", basename(untarred_16S)))
      if(length(untarred_ITS) > 0) {
        renamed <- file.rename(untarred_ITS, untarred_ITS_rename_to)
        if(!all(renamed)) {
          warning(sum(!renamed), " untarred ITS file(s) were not successfully renamed in sequencer run ", runID)
        }
        files_organized <- c(files_organized, basename(untarred_ITS[renamed]))
      }
      if(length(untarred_16S) > 0) {
        renamed <- file.rename(untarred_16S, untarred_16S_rename_to)
        if(!all(renamed)) {
          warning(sum(!renamed), " untarred 16S file(s) were not successfully renamed in sequencer run ", runID)
        }
        files_organized <- c(files_organized, basename(untarred_16S[renamed]))
      }

    # If the file is not tarred (typical case)
    } else {

      # rename file by appending sequencer run ID and moving to target gene-specific subdirectory
      if(grepl("_ITS[^/]*fastq", fn[i])) {
        rename_to <- file.path(outdir_sequence, "ITS", "0_raw", paste0("run", runID, "_", basename(fn[i])))
      } else if (grepl("_16S[^/]*fastq", fn[i])) {
        rename_to <- file.path(outdir_sequence, "16S", "0_raw", paste0("run", runID, "_", basename(fn[i])))
      } else {
        warning("Could not distinguish target gene of ", fn[i], ". We use regex patterns '_ITS[^/]*fastq' and '_16S[^/]*fastq'
              on the filename to make this distinction.")
        next
      }

      renamed <- file.rename(fn[i], rename_to)
      if(!renamed) {
        warning(fn[i], " was not successfully renamed.")
      } else {
        files_organized <- c(files_organized, basename(fn[i]))
      }
    }
  }
  if(length(files_organized) > 0) {
    message("Reorganized files can be found in ", file.path(outdir_sequence, "ITS"), " and ", file.path(outdir_sequence, "16S"), ".")
    return(files_organized)
  } else {
    warning("No files were successfully reorganized")
  }
}


#' Download NEON Soil Data Associated with Marker Gene Sequencing Data
#'
#' Downloads the following data products to the specified filepath:
#' (1) DP1.10078: "Soil chemical properties (Distributed periodic)";
#' (2) DP1.10086: "Soil physical properties (Distributed periodic)".
#' This function uses \code{\link[neonUtilities]{zipsByProduct}} to conduct the downloads.
#'
#' @param sites Either the string 'all', meaning all available sites, or a character vector of 4-letter NEON site codes, e.g. c('ONAQ','RMNP'). Defaults to PRESET_SITES parameter in params.R.
#' @param startYrMo,endYrMo Either NA, meaning all available dates, or a character vector in the form YYYY-MM, e.g. 2017-01. Defaults to PRESET_START_YR_MO in params.R.
#' @param outdir Location where output files are saved. Defaults to file.path(PRESET_OUTDIR, PRESET_OUTDIR_SOIL) in params.R.
#' @param checkFileSize TRUE or FALSE. Whether the user should be told the total file size before downloading. Defaults to PRESET_CHECK_FILE_SIZE in params.R.
#' @param return_data Whether to return metadata. If FALSE, simply downloads data to outdir and returns no value. Defaults to PRESET_RETURN_DATA.
#' @param overwrite TRUE or FALSE. If there is a previous download of this metadata in outdir, whether to overwrite that download.
#'
#' @return If return_data==TRUE, returns a dataframe consisting of joined soil data records from DP1.10078 ("Soil chemical properties (Distributed periodic)") and DP1.10086 ("Soil physical properties (Distributed periodic)"). Otherwise, no value is returned.
downloadRawSoilData <- function(sites = PRESET_SITES, startYrMo = PRESET_START_YR_MO, endYrMo = PRESET_END_YR_MO,
                                outdir = file.path(PRESET_OUTDIR, PRESET_OUTDIR_SOIL), checkFileSize = PRESET_CHECK_FILE_SIZE,
                                return_data = PRESET_RETURN_DATA, overwrite = FALSE) {
  # TODO: Why does this function sometimes return a few warnings of the form:
  # 1: In UseMethod("depth") :
  # no applicable method for 'depth' applied to an object of class "NULL"

  library(dplyr)
  library(neonUtilities)

  PRNUM_chem <- "10078"
  PRNUM_phys <- "10086"
  dpID_chem <- paste("DP1", PRNUM_chem, "001", sep=".")
  dpID_phys <- paste("DP1", PRNUM_phys, "001", sep=".")
  stackDir_chem <- file.path(outdir, paste0("filesToStack",PRNUM_chem))
  stackDir_phys <- file.path(outdir, paste0("filesToStack",PRNUM_phys))

  # If there is NO previous download at stackDir_chem
  if(!dir.exists(stackDir_chem)) {
    neonUtilities::zipsByProduct(dpID_chem, site=sites, startdate=startYrMo,
                                 enddate=endYrMo, package="expanded", check.size=checkFileSize, savepath = outdir)

    stackByTable(stackDir_chem, folder = TRUE)

    # Write site_and_date_range.txt file to record parameters
    write_sample_subset_params(stackDir_chem, sites, startYrMo, endYrMo)

    # If there IS a previous download at stackDir_chem
  } else {
    # If overwrite==TRUE,
    if(overwrite) {
      # delete existing directory
      unlink(stackDir_chem, recursive=TRUE)
      # and re-download data
      neonUtilities::zipsByProduct(dpID_chem, site=sites, startdate=startYrMo,
                                   enddate=endYrMo, package="expanded", check.size=checkFileSize,
                                   savepath=outdir)
      stackByTable(stackDir_chem, folder = TRUE)
      write_sample_subset_params(stackDir_chem, sites, startYrMo, endYrMo)

      # If overwrite==FALSE (default)
    } else {
      warn_already_downloaded(PRNUM_chem, outdir)
    }
  }

  # If there is NO previous download at stackDir_phys
  if(!dir.exists(stackDir_phys)) {
    neonUtilities::zipsByProduct(dpID_phys, site=sites, startdate=startYrMo,
                                 enddate=endYrMo, package="expanded", check.size=checkFileSize, savepath = outdir)

    stackByTable(stackDir_phys, folder = TRUE)

    # Write site_and_date_range.txt file to record parameters
    write_sample_subset_params(stackDir_phys, sites, startYrMo, endYrMo)

    # If there IS a previous download at stackDir_phys
  } else {
    # If overwrite==TRUE,
    if(overwrite) {
      # delete existing directory
      unlink(stackDir_phys, recursive=TRUE)
      # and re-download data
      neonUtilities::zipsByProduct(dpID_phys, site=sites, startdate=startYrMo,
                                   enddate=endYrMo, package="expanded", check.size=checkFileSize,
                                   savepath=outdir)
      stackByTable(stackDir_phys, folder = TRUE)
      write_sample_subset_params(stackDir_phys, sites, startYrMo, endYrMo)

      # If overwrite==FALSE (default)
    } else {
      warn_already_downloaded(PRNUM_phys, outdir)
    }
  }

  if(return_data) { # If FALSE, then simply downloads the data to the outdir

    path_soilCoreCollection <- file.path(stackDir_phys, "stackedFiles", "sls_soilCoreCollection.csv")
    path_soilMoisture <- file.path(stackDir_phys, "stackedFiles", "sls_soilMoisture.csv")
    path_soilpH <- file.path(stackDir_phys, "stackedFiles", "sls_soilpH.csv")
    # path_bgcSubsampling <- file.path(stackDir_phys, "stackedFiles", "sls_bgcSubsampling.csv")
    path_soilChemistry <- file.path(stackDir_chem, "stackedFiles", "sls_soilChemistry.csv")

    dat_soilCoreCollection <- read.delim(path_soilCoreCollection, sep=",", stringsAsFactors=FALSE)
    dat_soilMoisture <- read.delim(path_soilMoisture, sep=",", stringsAsFactors=FALSE)
    dat_soilpH <- read.delim(path_soilpH, sep=",", stringsAsFactors=FALSE)
    dat_soilChemistry <- read.delim(path_soilChemistry, sep=",", stringsAsFactors=FALSE)

    joining_cols <- c("domainID", "siteID", "plotID", "sampleID")

    # TODO: Confirm that the selected columns are sufficient for downstream analysis
    select(dat_soilCoreCollection, domainID, siteID, plotID, namedLocation, plotType, nlcdClass, coreCoordinateX, coreCoordinateY, geodeticDatum, decimalLatitude, decimalLongitude, elevation,
           sccSamplingProtocolVersion=samplingProtocolVersion, collectDate, sampleTiming, standingWaterDepth, nTransBoutType, sampleID, horizon, soilTemp, litterDepth, sampleTopDepth, sampleBottomDepth,
           soilSamplingDevice, geneticSampleID, sccDataQF=dataQF) %>%
      full_join(select(dat_soilMoisture, all_of(joining_cols), moistureSampleID, smSamplingProtocolVersion=samplingProtocolVersion, soilMoisture, smDataQF), by=joining_cols) %>%
      full_join(select(dat_soilpH, all_of(joining_cols), pHSampleID, pHSamplingProtocolVersion=samplingProtocolVersion, soilInWaterpH, soilInCaClpH, pHDataQF), by=joining_cols) %>%
      full_join(select(dat_soilChemistry, all_of(joining_cols), cnSampleID, nitrogenPercent, organicCPercent, CNratio, cnTestMethod=testMethod, cnInstrument=instrument, cnDataQF=dataQF), by=joining_cols) ->
      dat_soil
    # TODO: Filter by nTransBoutType to remove incubated samples?

    return(dat_soil)
  }
}


#' Download Sequence Metadata
#'
#' Loads soil marker gene sequencing metadata for specified target gene, site(s) and date(s),
#' with an option to download output by providing a valid output directory.
#'
#' Function by Lee F. Stanish and Clara Qin (2020). Currently available for testing only.
#'
#' @param sites Either the string 'all', meaning all available sites, or a character vector of 4-letter NEON site codes, e.g. c('ONAQ','RMNP'). Defaults to PRESET_SITES parameter in params.R.
#' @param startYrMo,endYrMo Either NA, meaning all available dates, or a character vector in the form YYYY-MM, e.g. 2017-01.
#' @param targetGene '16S' or 'ITS'.
#' @param sequencingRuns Either the string 'all', meaning all available sequencing runs, or a character vector of NEON sequencing run IDs, e.g. c('C25G9', 'B69PP').
#' @param dpID NEON data product of interest. Default is soil marker gene sequences, and currently code only works for this dpID.
#' @param outDir (Optional) If a local copy of the filtered metadata is desired, provide path to output directory.
#'
#' @return Data frame containing joined records from across the NEON soil marker gene sequence metadata, subsetted according to function arguments.
#' @examples
#' \dontrun{
#' meta <- downloadSequenceMetadataRev('all', '2015-01', '2016-01', '16S') # metadata is not saved to local directory
#' meta <- downloadSequenceMetadataRev('all', '2015-01', '2016-01', '16S', dir='./data/') # metadata is saved to local directory
#' }
downloadSequenceMetadataRev <- function(sites='all', startYrMo, endYrMo, targetGene= "all",
                                        sequencingRuns = "", dpID = "DP1.10108.001", outDir="") {
  # author: Lee Stanish
  # date: 2020-08-13
  # function loads soil marker gene sequencing metadata for target gene, site(s) and date(s)
  # option to download output by providing a valid output directory
  # sites: character vector of valid site ID's, or 'all' for all sites
  # targetGene: '16S',  'ITS', 'all'
  # startYrMo: start date, format YYYY-MM
  # endYrMo: end date, format YYYY-MM
  # dpID: NEON data product of interest. Default is soil marker gene sequences, and currently code only works for this dpID
  # outDir (optional): If a local copy of the filtered metadata is desired, provide path to output dir

  library(neonUtilities)
  library(plyr)
  library(dplyr)

  # check valid data values entered
  ## validate dpID ##
  if(!grepl("DP1", dpID) | !grepl('\\.001', dpID) | !grepl('10108|20280|20282', dpID)) {
    message("Invalid Data Product ID: must follow convention 'DP1.[5-digit value].001' and must be a marker genes data product ID")
    return(NULL)
  } else {
    dpID <- dpID
  }

  # validate target gene
  if(!grepl("16S|ITS|all", targetGene)) {
    message("Invalid targetGene: must be either '16S', 'ITS', 'all' ")
    return(NULL)
  } else {
    targetGene <- targetGene
  }

  # validate site(s)
  terrSiteList <- c("all","HARV","SCBI","OSBS","GUAN","UNDE","KONZ","ORNL","TALL","WOOD","CPER","CLBJ","YELL","NIWO",
                    "SRER","ONAQ","WREF","SJER","TOOL","BONA","PUUM","BART","BLAN","SERC","SCBI","DSNY","JERC","LAJA",
                    "TREE","STEI","KONA","UKFS","MLBS","GRSM","LENO","DELA","NOGP","DCFS","STER","RMNP","OAES","MOAB",
                    "JORN","ABBY","TEAK","SOAP","BARR","DEJU","HEAL")
  if(!any(sites %in% terrSiteList)){
    message("Invalid site(s): must be a valid NEON site or 'all'")
    return(NULL)
  } else {
    sites <- sites
  }

  message("loading metadata...")
  mmgL1 <- loadByProduct(dpID, sites, package = 'expanded', check.size = F, startdate = startYrMo, enddate = endYrMo) # output is a list of each metadata file


  # for target data product and targetGene: extract lists into data.frames
  if(grepl("10108", dpID)) {
    seq16S <- mmgL1$mmg_soilMarkerGeneSequencing_16S
    seqITS <- mmgL1$mmg_soilMarkerGeneSequencing_ITS
    raw <- mmgL1$mmg_soilRawDataFiles
    dna <- mmgL1$mmg_soilDnaExtraction
    seq <- rbind(seq16S, seqITS)

    if(targetGene=="16S") {
      message("filtering to 16S data")
      seq <- seq16S
    }
    if(targetGene=="ITS") {
      message("filtering to ITS data")
      seq <- seqITS
      raw <- mmgL1$mmg_soilRawDataFiles
    }
  }

  if(grepl("20280", dpID)) {
    seq16S <- mmgL1$mmg_benthicMarkerGeneSequencing_16S
    seqITS <- mmgL1$mmg_benthicMarkerGeneSequencing_ITS
    raw <- mmgL1$mmg_benthicRawDataFiles
    dna <- mmgL1$mmg_benthicDnaExtraction
    seq <- rbind(seq16S, seqITS)

    if(targetGene=="16S") {
      message("filtering to 16S data")
      seq <- seq16S
    }
    if(targetGene=="ITS") {
      message("filtering to ITS data")
      seq <- seqITS
    }
  }

  if(grepl("20282", dpID)) {
    seq16S <- mmgL1$mmg_swMarkerGeneSequencing_16S
    seqITS <- mmgL1$mmg_swMarkerGeneSequencing_ITS
    raw <- mmgL1$mmg_swRawDataFiles
    seq <- rbind(seq16S, seqITS)
    dna <- mmgL1$mmg_swDnaExtraction

    if(targetGene=="16S") {
      message("filtering to 16S data")
      seq <- seq16S
    }
    if(targetGene=="ITS") {
      message("filtering to ITS data")
      seq <- seqITS
    }
  }

  # remove unnecessary/redundant columns from tables
  raw <- select(raw, -domainID, -siteID, -namedLocation, -laboratoryName, -sequencingFacilityID, -collectDate, -dnaSampleCode)
  dna <- select(dna, -domainID, -siteID, -namedLocation, -laboratoryName, -collectDate)

  # convert factors to characters (bug in output of loadByProduct)
  i <- sapply(seq, is.factor)
  seq[i] <- lapply(seq[i], as.character)
  j <- sapply(raw, is.factor)
  raw[j] <- lapply(raw[j], as.character)
  j <- sapply(dna, is.factor)
  raw[j] <- lapply(dna[j], as.character)


  # If specified, filter by sequencing run ID
  if(sequencingRuns[1] != "") {
    raw <- raw[which(raw$sequencerRunID %in% sequencingRuns), ]
    # Validate sequencing run ID argument
    if(nrow(raw) == 0) {
      stop("After filtering by specified sequencing run ID(s), no records remain. Double-check your sequencing run  ID(s).")
    }
  }

  # Join sequencing metadata with raw data files metadata
  if(targetGene=="16S") {
    if(any(grepl("ITS", raw$rawDataFileName))) {
      rawCleaned <- raw[-grep("ITS", raw$rawDataFileName), ]
    } else {
      rawCleaned <- raw
    }
    joinedTarget <- left_join(rawCleaned, seq, by=c('dnaSampleID', 'sequencerRunID', 'internalLabID'))
    out <- joinedTarget[!is.na(joinedTarget$uid.y), ]
  }
  if(targetGene=="ITS") {
    if(any(grepl("16S", raw$rawDataFileName))) {
      rawCleaned <- raw[-grep("16S", raw$rawDataFileName), ]
    } else {
      rawCleaned <- raw
    }
    joinedTarget <- left_join(rawCleaned, seq, by=c('dnaSampleID', 'sequencerRunID', 'internalLabID'))
    out <- joinedTarget[!is.na(joinedTarget$uid.y), ]
  }
  if(targetGene=="all") {
    joinedTarget <- left_join(raw, seq, by=c('dnaSampleID', 'sequencerRunID', 'internalLabID'))
    out <- joinedTarget[!is.na(joinedTarget$uid.y), ]
    message(paste0(length(grep("16S", joinedTarget$rawDataFileName)), " 16S records and ", length(grep("ITS", joinedTarget$rawDataFileName)), " ITS records found."))
  }

  # clean up redundant column names
  names(out) <- gsub("\\.x", ".rawFiles", names(out))
  names(out) <- gsub("\\.y", ".seq", names(out))

  # join with DNA extraction metadata
  outDNA <- left_join(out, dna, by=c('plotID', 'dnaSampleID', 'internalLabID'))
  # clean up redundant column names
  names(outDNA) <- gsub("\\.x", ".seq", names(outDNA))
  names(outDNA) <- gsub("\\.y", ".dna", names(outDNA))
  names(outDNA)[names(outDNA)=="uid"] <- 'uid.dna'
  names(outDNA)[names(outDNA)=="remarks"] <- 'remarks.dna'
  names(outDNA)[names(outDNA)=="dataQF"] <- 'dataQF.dna'

  # download local copy if user provided output dir path
  if(outDir != "") {
    if(!dir.exists(outDir)) {
      dir.create(outDir)
    }
    write.csv(outDNA, paste0(outDir, "/mmg_soilMetadata_", targetGene, "_", Sys.Date(), ".csv"),
              row.names=F)
    print(paste0("metadata downloaded to: ", outDir, "/mmg_soilMetadata_", targetGene, "_", Sys.Date(), ".csv") )
  }
  return(outDNA)
}


#' Trim Primers from 16S Sequences
#'
#' Trims primers from 16S sequences using \code{\link[dada2]{filterAndTrim}}. This function
#' assumes that each read begins with its full primer sequence and operates
#' by truncating the beginning of each read by the length of its primer
#' sequence.
#'
#' @param fnFs Full names of fastq files containing forward-read sequence data.
#' @param fnRs Full names of fastq files containing reverse-read sequence data.
#' @param PATH_CUT Output directory. If it does not exist, it will be created.
#' @param PRIMER_16S_FWD DNA sequence of forward-read primer.
#' @param PRIMER_16S_REV DNA sequence of reverse-read primer.
#' @param MULTITHREAD Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#'
#' @return Integer matrix denoting the number of reads remaining after primer-trimming for each input file.
#'
#' @examples
#' \dontrun{
#' trimPrimers16S(c("sample1_R1.fastq", "sample2_R1.fastq"), c("sample1_R2.fastq", "sample2_R2.fastq"), "path/to/output", "CCTACGGGNBGCASCAG", "GACTACNVGGGTATCTAATCC", MULTITHREAD = TRUE)
#' }
trimPrimers16S <- function(fnFs, fnRs, PATH_CUT, PRIMER_16S_FWD, PRIMER_16S_REV, MULTITHREAD#, quiet=T
){

  # Create output directory
  if(!dir.exists(PATH_CUT)) dir.create(PATH_CUT)

  # Create primer-trimmed file paths in output directory
  fnFs.cut <- file.path(PATH_CUT, basename(fnFs))
  fnRs.cut <- file.path(PATH_CUT, basename(fnRs))

  # Just trim the length of fwd and rev primers
  qa.out <- filterAndTrim(fnFs, fnFs.cut, fnRs, fnRs.cut, multithread = MULTITHREAD,
                          matchIDs = TRUE, trimLeft = c(nchar(PRIMER_16S_FWD), nchar(PRIMER_16S_REV)), compress=F)

  colnames(qa.out) <- c("input", "trimmed")
  return(qa.out)
}


#' Trim Primers from ITS Sequences
#'
#' Trims primers from ITS sequences using cutadapt. Cutadapt must be installed in order for this to work. Currently only supports R1 (forward-read) files.
#'
#' @param fn Full names of input fastq files, including directory paths. Files that do not exist will be ignored; however, if all files do not exist, this function will throw an error. It is assumed that these are R1 (forward-read) files only.
#' @param path_cut Output directory. If it does not exist, it will be created.
#' @param primer_ITS_fwd,primer_ITS_rev Default PRIMER_ITS_FWD and PRIMER_ITS_REV in params.R. DNA sequence of forward-read primer and reverse-read primer, respectively.
#' @param cutadapt_path Default CUTADAPT_PATH in params.R. Path to cutadapt on your file system.
#' @param post_samplename_pattern1,post_samplename_pattern2 (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#' @param very_verbose Default FALSE. Whether to print output from cutadapt. Unlike some other "verbose" arguments associated with the functions in this package, this does not default to VERBOSE in params.R.
#' @param discard_untrimmed Default FALSE. Whether to discard reads where a primer could not be found, leaving only those reads in which a primer has been trimmed.
#'
#' @return No value is returned.
#'
#' @examples
#' \dontrun{
#' trimPrimersITS(c("sample1_ITS_R1.fastq", "sample2_ITS_R1.fastq"), "path/to/output", "CTTGGTCATTTAGAGGAAGTAA", multithread = TRUE)
#' }
trimPrimersITS <- function(fn, path_cut, primer_ITS_fwd = PRIMER_ITS_FWD, primer_ITS_rev = PRIMER_ITS_REV, cutadapt_path = CUTADAPT_PATH, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq", very_verbose=FALSE, discard_untrimmed=FALSE#, quiet=T
){

  fnFs <- fn[file.exists(fn) & grepl(post_samplename_pattern1, fn)]
  if(length(fnFs) == 0) stop("No files found at specified location. Check file path.")

  # Create output directory
  if(!dir.exists(path_cut)) dir.create(path_cut)

  # Create primer-trimmed file paths in output directory
  fnFs.cut <- file.path(path_cut, basename(fnFs))

  # Get reverse-complement sequence of reverse primer
  rev_rc <- dada2:::rc(primer_ITS_rev)

  # Trim forward primer and the reverse-complement of reverse primer off of R1 (forward reads)
  R1.flags <- paste("-g", primer_ITS_fwd, "-a", rev_rc)

  # Discard untrimmed reads?
  if(discard_untrimmed) {
    discard_untrimmed_flag <- "--discard-untrimmed"
  } else {
    discard_untrimmed_flag <- ""
  }

  # Run Cutadapt
  for(i in seq_along(fnFs)) {
    system2(cutadapt_path, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                                    "-o", fnFs.cut[i], # output files
                                    fnFs[i], # input files
                                    "--minimum-length", "1",
                                    discard_untrimmed_flag), # min length of cutadapted reads: >0
            stdout = ifelse(very_verbose, "", FALSE))
  }
}



#' Remove Unmatched Fastq Files
#'
#' Removes any forward-read files that do not have reverse-read counterparts, and vise versa.
#' This function is necessary because dada2::filterAndTrim() will throw an error
#' if the forward-read files and the reverse-read files are mismatched.
#'
#' @param fnFs Full name(s) of fastq file(s) containing forward-read sequence data.
#' @param fnRs Full name(s) of fastq file(s) containing reverse-read sequence data.
#' @param post_samplename_pattern (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#'
#' @return List of length 2. The first element is a vector of forward-read files that have reverse-read counterparts; the second element is a vector of reverse-read files that have forward-read counterparts.
#'
#' @examples
#' matched_fn <- remove_unmatched_files(c("sample1_R1.fastq", "sample2_R1.fastq"), c("sample1_R2.fastq", "sample2_R2.fastq", "sample3_R2.fastq"))
#' fnFs <- matched_fn[[1]]
#' fnRs <- matched_fn[[2]]
remove_unmatched_files <- function(fnFs, fnRs, post_samplename_pattern = "_R(1|2).*\\.fastq"){
  basefilenames_Fs <- sapply(strsplit(fnFs, post_samplename_pattern), `[`, 1)
  basefilenames_Rs <- sapply(strsplit(fnRs, post_samplename_pattern), `[`, 1)
  rm_from_fnFs <- which(!(basefilenames_Fs %in% basefilenames_Rs))
  rm_from_fnRs <- which(!(basefilenames_Rs %in% basefilenames_Fs))

  if (length(c(rm_from_fnFs, rm_from_fnRs)) == 0){
    if(VERBOSE == TRUE) message("All R1/R2 files had counterparts.\n")
  } else {
    for(i in rm_from_fnFs) {
      if(VERBOSE) message(paste(basefilenames_Fs[i], "does not have an R2 counterpart. Omitting from this analysis.\n"))
      fnFs <- fnFs[-rm_from_fnFs]
    }
    for(i in rm_from_fnRs) {
      if(VERBOSE) message(paste(basefilenames_Rs[i], "does not have an R1 counterpart. Omitting from this analysis.\n"))
      fnRs <- fnRs[-rm_from_fnRs]
    }
  }
  return(list(R1=fnFs, R2=fnRs))
}


#' Filter 16S Sequences
#'
#' Applies a quality filter to 16S sequence fastq files via the \code{\link[dada2]{filterAndTrim}} function.
#' It is assumed that both forward- and reverse-read files are included.
#'
#' @param PATH_CUT Path to directory containing input fastq files.
#' @param PATH_FILTERED Path to output directory where filtered fastq files will be written.
#' @param MULTITHREAD Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#' @param MAX_EE_FWD,MAX_EE_REV The maximum number of expected errors allowed in forward/reverse reads. Read-pairs with expected errors that exceed this threshold will be removed.
#' @param TRUNC.LENGTHS Default NULL. Single integer: truncation length to use across all files. Two-integer vector: truncation length to use for the forward-read and reverse-read files, respectively. If NULL (default), determines truncation length(s) based on \code{\link{getTruncationLength}} with a quality score threshold of trunc_qscore.
#' @param trunc_qscore Default 23. Quality score at which point to truncate each read, if TRUNC.LENGTHS is null.
#' @param post_samplename_pattern1,post_samplename_pattern2 (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#'
#' @return Two-column matrix displaying the number of reads in input vs. output for each file.
qualityFilter16S <- function(PATH_CUT, PATH_FILTERED, MULTITHREAD, MAX_EE_FWD, MAX_EE_REV, TRUNC.LENGTHS = NULL, trunc_qscore = 23, post_samplename_pattern = "_R(1|2).fastq" #, mean = FALSE
){
  fnFs <- sort(list.files(PATH_CUT, pattern = "_R1.fastq", full.names = TRUE)) # primer-trimmed
  fnRs <- sort(list.files(PATH_CUT, pattern = "_R2.fastq", full.names = TRUE)) # primer-trimmed

  sample.names <- sapply(strsplit(basename(fnFs), post_samplename_pattern), `[`, 1) # extract sample names
  filtFs <- file.path(PATH_FILTERED, paste0(sample.names, "_filt_R1.fastq.gz")) # create filtered filenames
  filtRs <- file.path(PATH_FILTERED, paste0(sample.names, "_filt_R2.fastq.gz")) # create filtered filenames

  # CLARA: I think it will be easier to keep track of samples across processing steps if we
  # keep the filenames the same, and only vary the directory location at each step.

  if (is.null(TRUNC.LENGTHS)){
    # GETTING TRUNCATION LENGTH USING A SUBSET OF QUALITY SCORES
    n_files <- min(length(unique(sample.names)), 30)

    fwd.trunc.lengths <- list()
    rev.trunc.lengths <- list()
    fwd.trunc.lengths[1:n_files] <- getTruncationLength(fnFs[1:n_files], verbose = VERBOSE,
                                                        qscore = trunc_qscore)
    rev.trunc.lengths[1:n_files] <- getTruncationLength(fnRs[1:n_files], verbose = VERBOSE,
                                                        qscore = trunc_qscore)
    # TODO: Why this specific quality score?

    # remove any samples that have low-quality reads early on, to avoid spoiling the whole run
    # idk how to loop this command
    to_remove <- which(unlist(fwd.trunc.lengths) < 100 | unlist(rev.trunc.lengths) < 100)
    if (length(to_remove) > 0){
      fwd.trunc.lengths <- fwd.trunc.lengths[-to_remove]
      rev.trunc.lengths <- rev.trunc.lengths[-to_remove]
    }

    # if (mean == T){
    #   fwd.trunc.length <- round(mean(unlist(fwd.trunc.lengths))) # tested this with mean vs minimum quality-score threshold;
    #   rev.trunc.length <- round(mean(unlist(rev.trunc.lengths))) # minimum retained more reads after filtering.
    # } else {
    fwd.trunc.length <- round(min(unlist(fwd.trunc.lengths))) # tested this with mean vs minimum quality-score threshold;
    rev.trunc.length <- round(min(unlist(rev.trunc.lengths))) # minimum retained more reads after filtering.
    # }

    # Set the minimum lengths for fwd/reverse reads (if they're too short, they cannot be merged).
    # TODO: What was the reasoning for selecting these specific lengths?
    if (rev.trunc.length < 200) rev.trunc.length <- 200
    if (fwd.trunc.length < 245) fwd.trunc.length <- 245
    cat(paste0("Fwd truncation length: ", fwd.trunc.length, "\nRev truncation length: ", rev.trunc.length, "\n"))
    TRUNC.LENGTHS <- c(fwd.trunc.length, rev.trunc.length)
  }

  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                       multithread=MULTITHREAD,
                       truncLen=TRUNC.LENGTHS,
                       maxEE = c(MAX_EE_FWD, MAX_EE_REV),
                       truncQ=TRUNC_Q, matchIDs = TRUE,
                       compress=TRUE, maxN=0)

  return(out)
}


# TODO: How can we add "..." to pass arbitrary parameters to filterAndTrim?

#' Filter ITS Sequences
#'
#' Applies a quality filter to ITS sequence fastq files via the \code{\link[dada2]{filterAndTrim}} function.
#' Currently only supports filtering forward-read sequences.
#'
#' @param fn Full names of input fastq files, including directory paths. Files that do not exist will be ignored; however, if all files do not exist, this function will throw an error. It is assumed that these are R1 (forward-read) files only.
#' @param path_filtered Path to output directory where filtered fastq files will be written.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#' @param maxEE Default Inf. The maximum number of expected errors allowed in forward reads. Read-pairs with expected errors that exceed this threshold will be removed.
#' @param truncQ Default 2. Quality score at which point to truncate each read.
#' @param minLen Default 20. Read with length less than minLen are removed. minLen is enforced after trimming and truncation.
#' @param maxN Default 0. Maximum number of ambiguous bases "N" to allow in forward reads.
#' @param post_samplename_pattern1,post_samplename_pattern2 (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#' @param trimLeft (Optional). Default 0. The number of nucleotides to remove from the start of each read.
#'
#' @return Two-column matrix displaying the number of reads in input vs. output for each file.
qualityFilterITS <- function(fn, path_filtered, multithread = MULTITHREAD, maxEE = Inf, truncQ = 2, minLen = 20, maxN = 0, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq", trimLeft = 0 #, mean = FALSE
){
  fnFs <- fn[file.exists(fn) & grepl(post_samplename_pattern1, fn)]
  if(length(fnFs) == 0) stop("No files found at specified location. Check file path.")

  sample.names <- sapply(strsplit(basename(fnFs), paste(post_samplename_pattern1, "|", post_samplename_pattern2)), `[`, 1) # extract sample names
  # filtFs <- file.path(path_filtered, paste0(sample.names, "_filt_R1.fastq.gz")) # create filtered filenames
  # filtRs <- file.path(path_filtered, paste0(sample.names, "_filt_R2.fastq.gz")) # create filtered filenames

  filtFs <- file.path(path_filtered, basename(fnFs))

  out <- filterAndTrim(fnFs, filtFs,
                       multithread = multithread,
                       maxEE = maxEE,
                       truncQ = truncQ,
                       minLen = minLen,
                       compress = TRUE,
                       maxN = maxN,
                       trimLeft = trimLeft)

  rownames(out) <- sample.names

  return(out)
}

#' Get Truncation Length
#'
#' Decides on truncation length for trimmed reads based on quality score means. Default cutoff is a score of 30. Warns you if the beginning (first 10) bases are low-quality, but returns the first low-quality base after the 10th base. If no bases are below score, returns the last base. The truncation lengths can be aggregated (e.g. minimum) and used as an argument in \code{\link[dada2]{filterAndTrim}}.
#'
#' @param fl Full names of fastq files.
#' @param qscore Default 30. Mean quality score threshold at which to truncate the remainder of the read.
#' @param n Default 5e+05. The number of reads to sample when processing fastq files.
#' @param verbose Default TRUE. Whether to return message regarding truncation length for each file. Includes warning messages.
#'
#' @return Integer vector of truncation lengths to use for the set of fastq files.
#'
#' @examples
#' \dontrun{
#' trunc_len_fwd <- getTruncationLength(c("sample1_R1.fastq", "sample2_R1.fastq"))
#' trunc_len_rev <- getTruncationLength(c("sample1_R2.fastq", "sample2_R2.fastq"))
#' }
getTruncationLength <-function (fl, qscore = 30, n = 5e+05, verbose = TRUE){
  trunc_lengths <- data.frame(file = character(0), early_lowqual = numeric(0),
                              trunc_length = numeric(0))
  for (f in fl) {
    srqa <- ShortRead::qa(f, n = n)
    df <- srqa[["perCycle"]]$quality
    means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count,
                                                          df$Cycle)
    lowqual <- which(means < qscore)
    n_lowqual_in_first_20 <- sum(lowqual <= 20)
    first_lowqual_after_20 <- lowqual[lowqual > 20][1]

    trunc_lengths <- rbind(trunc_lengths, data.frame(file = f,
                                                     early_lowqual = n_lowqual_in_first_20,
                                                     trunc_length = ifelse(is.na(first_lowqual_after_20), length(means), first_lowqual_after_20)))

    if (verbose == TRUE) {
      if(lowqual_in_first_20 > 0){
        cat(paste0("Note: for sample ", basename(f),': ', n_lowqual_in_first_20, ' base(s) within first 20 bases are below your quality score. Consider trimming left.\n'))
      }
      if (is.na(lowqual_after_20)){
        cat(paste0(basename(f), ': After first 20 bases, no bases with a mean under your quality score found. Truncate at end of read, base: ', length(means),'\n'))

      } else if (!is.na(lowqual_after_20)){
        cat(paste0(basename(f),': After first 20 bases, first mean below your quality score is base: ',first_lowqual_after_20,'\n'))

      } else "Something's wrong. Inspect function/reads."
    } # end printing to console
  } # end loop
  return(trunc_lengths$trunc_length)
}



#' Run Dada on Paired-End 16S Sequences
#'
#' Runs the Dada algorithm to infer sample composition from paired-end 16S fastq files.
#' This implementation is based on Ben Callahan's vignette at \url{https://benjjneb.github.io/dada2/bigdata_paired.html}.
#'
#' @param PATH_FILTERED Path to directory containing input fastq files.
#' @param MULTITHREAD Whether to use multithreading.
#' @param VERBOSE Default FALSE. Whether to print messages regarding the dimensions of the resulting sequence table and the distribution of sequence lengths.
#' @param seed (Optional) Integer to use as random seed for reproducibility.
#'
#' @return A list of three elements. \strong{seqtab} is the sequence table before removing chimeras, \strong{seqtab.nochim} is the sequence table after removing chimeras, and \strong{track} is a data frame displaying the number of reads remaining for each sample at various points throughout the processing pipeline.
#'
#' @examples
#' \dontrun{
#' seqtab.list <- runDada16S('./seq/filtered/', TRUE, TRUE, 1010100)
#' }
runDada16S <- function(PATH_FILTERED, MULTITHREAD, VERBOSE = FALSE, seed = NULL){
  if (!is.null(seed)) set.seed(seed)

  # File parsing
  filtFs <- list.files(PATH_FILTERED, pattern="_R1.fastq", full.names = TRUE)
  filtRs <- list.files(PATH_FILTERED, pattern="_R2.fastq", full.names = TRUE)

  # Keep files with counterpart
  remove_unmatched_files(filtFs, filtRs)

  # Create outnames
  sample.names <- sapply(strsplit(basename(filtFs), "_R(1|2).fastq"), `[`, 1) # Assumes filename = samplename_RX.fastq.gz
  sample.namesR <- sapply(strsplit(basename(filtRs), "_R(1|2).fastq"), `[`, 1) # Assumes filename = samplename_RX.fastq.gz
  if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  # Learn forward and reverse error rates
  errF <- learnErrors(filtFs, nbases=1e7, randomize=TRUE, multithread=MULTITHREAD)
  errR <- learnErrors(filtRs, nbases=1e7, randomize=TRUE, multithread=MULTITHREAD)

  # Create output vectors
  mergers <- vector("list", length(sample.names))
  derepF <- vector("list", length(sample.names))
  derepR <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  names(derepF) <- sample.names
  names(derepR) <- sample.names

  # Sample inference and merger of paired-end reads
  for(i in 1:length(sample.names)) {
    sam <- sample.names[[i]]
    cat("Processing:", sam, "\n")

    derepF[[i]] <- derepFastq(filtFs[[i]])
    derepR[[i]] <- derepFastq(filtRs[[i]])

    ddF <- dada(derepF[[i]], err=errF, multithread=TRUE)
    ddR <- dada(derepR[[i]], err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF[[i]], ddR, derepR[[i]], maxMismatch=1, minOverlap = 6)
    mergers[[i]] <- merger
    cat(paste("Exact sequence variants inferred for sample:", sam,". \n"))
  }
  #rm(derepF); rm(derepR)
  # Construct sequence table and remove chimeras
  seqtab <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=MULTITHREAD, verbose=VERBOSE)

  if(VERBOSE){
    cat(paste0("\n\nDimensions of ESV table: ", dim(seqtab.nochim)[1], " samples, ", dim(seqtab.nochim)[2], " ESVs\n"))
    cat("\nDistribution of sequence lengths (should be bimodal for v3v4 region):")
    print(table(nchar(getSequences(seqtab))))
    cat(paste0("\n", round(sum(seqtab.nochim)/sum(seqtab), 3), " reads remain after removal of chimeras"))
  }
  if(nrow(seqtab.nochim) > 1) {
    track <- cbind.data.frame(derepF = sapply(derepF, getN),
                              derepR = sapply(derepR, getN),
                              denoisedF = sapply(ddF, getN),
                              denoisedR = sapply(ddR, getN),
                              merged = sapply(mergers, getN),
                              nonchim = rowSums(seqtab.nochim))
  # If processing a single sample, remove the sapply calls: e.g. replace
  # sapply(dadaFs, getN) with getN(dadaFs)
  } else {
    track <- cbind.data.frame(derepF = getN(derepF),
                              derepR = getN(derepR),
                              denoisedF = getN(ddF),
                              denoisedR = getN(ddR),
                              merged = getN(mergers),
                              nonchim = sum(seqtab.nochim))
  }
  return(list("seqtab" = seqtab, "seqtab.nochim" = seqtab.nochim, "track" = track))
}


#' Run Dada on R1 ITS sequences
#'
#' Runs the Dada algorithm to infer sample composition from R1 ITS fastq files.
#' We recommend using only the R1 reads because attempting to merge with R2
#' reads may result in less accurate representations of the fungal community
#' composition (Pauvert et al., 2019).
#' This implementation is based on Ben Callahan's vignettes at \url{https://benjjneb.github.io/dada2/bigdata.html}
#' and \url{https://benjjneb.github.io/dada2/ITS_workflow.html}.
#'
#' @param fn Full names of input fastq files, including directory paths. Files that do not exist will be ignored; however, if all files do not exist, this function will throw an error. It is assumed that these are R1 (forward-read) files only, and that they contain no primers, and that they have been filtered.
#' @param multithread Whether to use multithreading.
#' @param verbose Default FALSE. Whether to print messages regarding the dereplication step, the denoising step, and the dimensions of the resulting sequence table and the distribution of sequence lengths.
#' @param seed (Optional) Integer to use as random seed for reproducibility.
#' @param post_samplename_pattern1,post_samplename_pattern2 (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#'
#' @return A list of three elements. \strong{seqtab} is the sequence table before removing chimeras, \strong{seqtab.nochim} is the sequence table after removing chimeras, and \strong{track} is a data frame displaying the number of reads remaining for each sample at various points throughout the processing pipeline.
#'
#' @examples
#' \dontrun{
#' seqtab.list <- runDadaITS('./seq/filtered/', TRUE, TRUE, 1010100)
#' }
runDadaITS <- function(fn, multithread, verbose = FALSE, seed = NULL, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq"){
  if (!is.null(seed)) set.seed(seed)

  # File parsing
  filtFs <- fn[file.exists(fn) & grepl(post_samplename_pattern1, fn)]
  if(length(filtFs) == 0) stop("No files found at specified location. Check file path.")

  # Create outnames
  sample.names <- sapply(strsplit(basename(filtFs), paste(post_samplename_pattern1, "|", post_samplename_pattern2)), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
  names(filtFs) <- sample.names

  # Learn error rates
  errF <- learnErrors(filtFs, nbases=1e7, randomize=TRUE, multithread=multithread)

  # Create output vector
  derepF <- vector("list", length(sample.names))
  ddF <- vector("list", length(sample.names))
  names(derepF) <- sample.names
  names(ddF) <- sample.names

  # Sample inference and merger of paired-end reads
  for(i in 1:length(sample.names)) {
    sam <- sample.names[[i]]
    cat("Processing:", sam, "\n")
    derepF[[i]] <- derepFastq(filtFs[[i]], verbose=verbose)
    ddF[[i]] <- dada(derepF[[i]], err=errF, multithread=MULTITHREAD, verbose=verbose)
    cat(paste("Exact sequence variants inferred for sample:", sam,". \n"))
  }
  # rm(derepF);

  # Construct sequence table and remove chimeras
  seqtab <- makeSequenceTable(ddF)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=multithread, verbose=verbose)

  if(verbose){
    cat(paste0("\n\nDimensions of ESV table: ", dim(seqtab.nochim)[1], " samples, ", dim(seqtab.nochim)[2], " ESVs\n"))
    cat("\nDistribution of sequence lengths:")
    print(table(nchar(getSequences(seqtab))))
    cat(paste0("\n", round(sum(seqtab.nochim)/sum(seqtab), 3), " reads remain after removal of chimeras"))
  }

  track <- cbind.data.frame(derepF = sapply(derepF, getN),
                            denoisedF = sapply(ddF, getN),
                            nonchim = rowSums(seqtab.nochim))
  return(list("seqtab" = seqtab, "seqtab.nochim" = seqtab.nochim, "track" = track))
}



#' Get Number of Reads
#'
#' Gets the total number of reads in a dada-class object.
#' Primarily used for tracking reads throughout the pipeline.
#'
#' @param x A dada-class object.
#'
#' @return Integer; total number of reads in the dada-class object.
getN <- function(x) {
  sum(getUniques(x))
}




####################################################################
# The following functions will be deprecated with the completion of
# downloadSequenceMetadataRev()

#' Write sample subsetting parameters
#'
#' Write site and date range parameters for subsetting samples. To be called within other functions.
#'
#' @param writeDir Directory where the subsetting parameters should be written.
#' @param sites Vector containing the names of sites in the subset.
#' @param startYrMo,endYrMo Dates describing the time range of subset, formatted as "YYYY-MM".
#' @param target_genes Marker gene used for sequencing. Currently must be one of "16S", "ITS".
#' @param sequencing_runs Vector naming the sequencing runs in the subset.
#'
#' @details This function writes the sample subsetting parameters to a file. In addition to
#' keeping a record, it is also meant to be used with \code{\link{warn_already_downloaded}}
#' to detect when a data product has already been downloaded, and if so, with what subsetting
#' parameters.
#'
#' @return No value is returned.
#'
#' @examples
#' \dontrun{
#' write_sample_subset_params(writeDir="./raw_data/", sites="all", startYrMo="2016-01", endYrMo="2019-01", target_genes="ITS", sequencing_runs=c("B69PP","C25G9"))
#' }
write_sample_subset_params <- function(writeDir, sites, startYrMo, endYrMo,
                                       target_genes="", sequencing_runs="") {
  # paste("writing in ", writeDir)
  write.table(
    data.frame(x1=c("sites", "startYrMo", "endYrMo", "target_genes", "sequencing_runs"),
               x2=c(paste0(sites, collapse=","), startYrMo, endYrMo, paste0(target_genes, collapse=","), paste0(sequencing_runs, collapse=","))),
    file = file.path(writeDir, "sample_subset_params.txt"),
    sep=":", col.names=FALSE, quote=FALSE, row.names=FALSE
  )
}

#' Warn That Metadata Is Already Downloaded
#'
#' Issue a warning that metadata has already been downloaded, and indicates the sites and date
#' ranges for which it was downloaded. To be called within other functions.
#'
#' @param PRNUM NEON data product ID, e.g. "10108".
#' @param outdir Directory where data product has been downloaded (not the stacked folder within it).
#'
#' @return No value is returned.
warn_already_downloaded <- function(PRNUM, outdir) {
  stackDir <- paste(outdir, paste0("filesToStack",PRNUM), sep="/")
  print(paste0("Warning: Data product ", PRNUM,
               " has already been downloaded to ", outdir, "."))
  if(file.exists(file.path(paste(stackDir, SAMPLE_SUBSET_PARAMS_FILENAME, sep="/")))) {
    site_and_date_range <- read.table(paste(stackDir, SAMPLE_SUBSET_PARAMS_FILENAME, sep="/"),
                                      header=FALSE, sep=":")
    print(paste0("Ensure that the following describes your intended data subset, ",
                 "or else stop the current process and re-run with overwrite==TRUE. ",
                 "  sites:", site_and_date_range[1,2],
                 "  startYrMo:", site_and_date_range[2,2],
                 "  endYrMo:", site_and_date_range[3,2],
                 "  target_genes:", site_and_date_range[4,2],
                 "  sequencing_runs:", site_and_date_range[5,2]))
  }
}


#' Download NEON Marker Gene Sequencing Metadata
#'
#' Downloads NEON data product DP1.10108 ("Soil microbe marker gene sequences"),
#' which contains information on how to access raw sequence data. This function
#' uses \code{\link[neonUtilities]{zipsByProduct}} to conduct the downloads.
#'
#' @param sites Either the string 'all', meaning all available sites, or a character vector of 4-letter NEON site codes, e.g. c('ONAQ','RMNP'). Defaults to PRESET_SITES parameter in params.R.
#' @param startYrMo,endYrMo Either NA, meaning all available dates, or a character vector in the form YYYY-MM, e.g. 2017-01. Defaults to PRESET_START_YR_MO in params.R.
#' @param outdir Location where output files are saved. Defaults to file.path(PRESET_OUTDIR, PRESET_OUTDIR_SEQMETA) in params.R.
#' @param checkFileSize TRUE or FALSE. Whether the user should be told the total file size before downloading. Defaults to PRESET_CHECK_FILE_SIZE in params.R.
#' @param return_data Whether to return part of the metadata: the table NEON table "mmg_soilRawDataFiles.csv" within NEON data product DP1.10108. If FALSE, simply downloads data to outdir and returns no value. Defaults to PRESET_RETURN_DATA.
#' @param target_genes 'ITS', '16S', or 'all'. Defaults to TARGET_GENE in params.R.
#' @param sequencing_runs Either the string 'all', meaning all available sequencing runs, or a character vector of NEON sequencing run IDs, e.g. c('C25G9', 'B69PP'). Defaults to SEQUENCING_RUNS parameter in params.R.
#' @param overwrite TRUE or FALSE. If there is a previous download of this metadata in outdir, whether to overwrite that download.
#'
#' @return If return_data==TRUE, returns the mmg_soilRawDataFiles table from NEON.DP1.10108 ("Soil microbe marker gene sequences"). Otherwise, no value is returned.
downloadSequenceMetadata <- function(sites = PRESET_SITES, startYrMo = PRESET_START_YR_MO, endYrMo = PRESET_END_YR_MO,
                                     outdir = file.path(PRESET_OUTDIR, PRESET_OUTDIR_SEQMETA), checkFileSize = PRESET_CHECK_FILE_SIZE,
                                     return_data = PRESET_RETURN_DATA, target_genes = TARGET_GENE, sequencing_runs = SEQUENCING_RUNS,
                                     overwrite = FALSE) {

  library(neonUtilities)

  PRNUM <- "10108"
  dpID <- paste("DP1", PRNUM, "001", sep=".")
  stackDir <- file.path(outdir, paste0("filesToStack",PRNUM))

  # If there is NO previous download at stackDir
  if(!dir.exists(stackDir)) {
    neonUtilities::zipsByProduct(dpID, site=sites, startdate=startYrMo,
                                 enddate=endYrMo, package="expanded", check.size=checkFileSize,
                                 savepath=outdir)

    # print(paste("Attempt to stackByTable in", stackDir))
    stackByTable(stackDir, folder = TRUE)

    # Write file that records parameters
    write_sample_subset_params(stackDir, sites, startYrMo, endYrMo, target_genes, sequencing_runs)

    # If there IS a previous download at stackDir
  } else {
    # If overwrite==TRUE,
    if(overwrite) {
      # delete existing directory
      unlink(stackDir, recursive=TRUE)
      # and re-download data
      neonUtilities::zipsByProduct(dpID, site=sites, startdate=startYrMo,
                                   enddate=endYrMo, package="expanded", check.size=checkFileSize,
                                   savepath=outdir)
      # print(paste("Attempt to stackByTable in", stackDir))
      stackByTable(stackDir, folder = TRUE)
      write_sample_subset_params(stackDir, sites, startYrMo, endYrMo, target_genes, sequencing_runs)

      # If overwrite==FALSE (default)
    } else {
      warn_already_downloaded(PRNUM, outdir)
    }
  }

  if(return_data) { # If not, then simply downloads the data to the outdir
    path <- paste(stackDir, "stackedFiles",
                  "mmg_soilRawDataFiles.csv", sep="/")
    dat <- read.delim(path, sep=",", stringsAsFactors=FALSE)

    return(dat)
  }
}
