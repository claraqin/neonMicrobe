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
                                  seqmeta_dirname=PRESET_OUTDIR_SEQMETA, soil_dirname=PRESET_OUTDIR_SOIL, 
																	taxref_dirname=PRESET_OUTDIR_TAXREF) {
  if(is.null(seq_dirname) | seq_dirname == "") seq_dirname <- "raw_sequence"
  if(is.null(seqmeta_dirname) | seqmeta_dirname == "") seqmeta_dirname <- "sequence_metadata"
  if(is.null(soil_dirname) | soil_dirname == "") soil_dirname <- "soil"

  seq_dir <- file.path(base_dir, seq_dirname)
  seqmeta_dir <- file.path(base_dir, seqmeta_dirname)
  soil_dir <- file.path(base_dir, soil_dirname)
  taxref_dir <- file.path(base_dir, taxref_dirname)
  
  message("Building output directories from base directory '", base_dir, "'.")
  message("Using '", seq_dir, "' as sequence subdirectory.")
  message("Using '", seqmeta_dir, "' as sequence metadata subdirectory.")
  message("Using '", soil_dir, "' as soil data subdirectory.")
  message("Using '", taxref_dir, "' as taxonomic reference data subdirectory.")
  
  seq_its_dir <- file.path(seq_dir, "ITS")
  seq_16s_dir <- file.path(seq_dir, "16S")

  # If preset output directories do not exist, create them
  if(!dir.exists(seq_dir)) dir.create(seq_dir, recursive=TRUE)
  if(!dir.exists(seqmeta_dir)) dir.create(seqmeta_dir, recursive=TRUE)
  if(!dir.exists(soil_dir)) dir.create(soil_dir, recursive=TRUE)
  if(!dir.exists(taxref_dir)) dir.create(taxref_dir, recursive=TRUE)
  
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
#' @return Returns (invisibly) a list of integer codes: 0 indicates success of downloads and a non-zero integer indicates failure. See the help page for \code{\link[utils]{download.file}} for more details.
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
      if(nrow(metadata) > 0) {
        message(length(tar_ind), " row(s) in metadata associated with batch-level sequence data
                were ignored prior to downloading raw sequence data.")
      } else {
        stop("No rows remain in metadata after removing rows associated with batch-level sequence
             data. Consider setting ignore_tar_files = FALSE")
      }
    }
  }

  u.urls <- unique(metadata$rawDataFilePath)
  fileNms <- gsub('^.*\\/', "", u.urls)
  print(paste("There are", length(u.urls), "unique raw sequence files to download."))

  download_success <- list()
  for(i in 1:length(u.urls)) {
    ## TODO: download.file commands should be replaced by neonUtilities::zipsByURI
    tryCatch({
      download_success[[i]] <- download.file(
        url = as.character(u.urls[i]),
        destfile = ifelse(dir.exists(outdir), paste(outdir, fileNms[i], sep="/"), paste(getwd(), fileNms[i], sep="/" )),
        quiet = !verbose)
    }, error = function(e) { # Occasionally an error arises because _fastq should be replaced by .fastq
      tryCatch({
        revised_url <- sub("_fastq", ".fastq", as.character(u.urls[i]))
        download_success[[i]] <- download.file(
          url = revised_url,
          destfile = ifelse(dir.exists(outdir), paste(outdir, fileNms[i], sep="/"), paste(getwd(), fileNms[i], sep="/" )),
          quiet = !verbose)
      }, error = function(f) {
        download_success[[i]] <- 2
      })
    })
    if(dir.exists(outdir)) {
      if(verbose) message("Finished downloading ", paste(outdir, fileNms[i], sep="/"))
    } else {
      if(verbose) message("Finished downloading ", paste(getwd(), fileNms[i], sep="/" ))
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
#' @param outdir_sequence Default file.path(PRESET_OUTDIR, PRESET_OUTDIR_SEQUENCE). Directory where raw sequence files can be found before reorganizing.
#' @param verbose If TRUE, prints message each time a file is reorganized.
#'
#' @return Character vector of the files (including files within tarballs) that were successfully reorganized. If no files were successfully reorganized, returns no value.
organizeRawSequenceData <- function(fn, metadata, outdir_sequence = file.path(PRESET_OUTDIR, PRESET_OUTDIR_SEQUENCE), verbose = TRUE) {
  library(R.utils)

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
      if(verbose) message("Reorganized ", basename(fn[i]))
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


#' Download NEON Soil Data Associated with Marker Gene Sequencing Data
#'
#' Downloads either or both of the following data products:
#' (1) DP1.10078.001: "Soil chemical properties (Distributed periodic)";
#' (2) DP1.10086.001: "Soil physical properties (Distributed periodic)".
#' This function uses \code{\link[neonUtilities]{loadByProduct}} to conduct the downloads.
#'
#' @param sites Either the string 'all', meaning all available sites, or a character vector of 4-letter NEON site codes, e.g. c('ONAQ','RMNP'). Defaults to PRESET_SITES parameter in params.R.
#' @param startYrMo,endYrMo Either NA, meaning all available dates, or a character vector in the form YYYY-MM, e.g. 2017-01. Defaults to PRESET_START_YR_MO in params.R.
#' @param dpID NEON data product(s) of interest. Default is both DP1.10078.001 ("Soil chemical properties (Distributed periodic)") and DP1.10086.001 ("Soil physical properties (Distributed periodic)").
#' @param outDir (Optional) If a local copy of the filtered metadata is desired, provide path to output directory.
#'
#' @return If return_data==TRUE, returns a dataframe consisting of joined soil data records from DP1.10078 ("Soil chemical properties (Distributed periodic)") and DP1.10086 ("Soil physical properties (Distributed periodic)"). Otherwise, no value is returned.
downloadRawSoilData <- function(sites='all', startYrMo, endYrMo,
                                dpID = c("DP1.10078.001", "DP1.10086.001"), outDir="") {

  library(dplyr)
  library(neonUtilities)

  # check valid data values entered
  ## validate dpID ##
  if(!all(grepl("DP1", dpID) & grepl('\\.001', dpID) & grepl('10078|10086', dpID))) {
    message("Invalid Data Product ID: must follow convention 'DP1.[5-digit value].001' and must be a distributed periodic soil data product ID")
    return(NULL)
  } else {
    dpID <- dpID
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

  slsL1 <- list()
  message("loading soil data...")
  for(i in 1:length(dpID)) {
    slsL1[[dpID[i]]] <- tryCatch({
      loadByProduct(dpID[i], sites, package = 'expanded', check.size = F, startdate = startYrMo, enddate = endYrMo) # output is a list of lists of each soil data file
    }, error = function(e) {
      warning("No data was found for data product ", dpID[i], " at the specified sites and dates.")
      NA
    })
  }

  joining_cols <- c("domainID", "siteID", "plotID", "sampleID")

  if(any(grepl('10086', dpID)) & !all(is.na(slsL1[["DP1.10086.001"]]))) {
    # start with soilCoreCollection data...
    dat_soil_phys <-
      dplyr::select(slsL1[["DP1.10086.001"]]$"sls_soilCoreCollection",
                    domainID, siteID, plotID, namedLocation, plotType, nlcdClass, coreCoordinateX, coreCoordinateY, geodeticDatum, decimalLatitude, decimalLongitude, elevation,
                    sccSamplingProtocolVersion=samplingProtocolVersion, collectDate, sampleTiming, standingWaterDepth, nTransBoutType, sampleID, horizon, soilTemp, litterDepth,
                    sampleTopDepth, sampleBottomDepth, soilSamplingDevice, geneticSampleID, sccDataQF=dataQF) %>%
      # merge with soilMoisture data...
      full_join(dplyr::select(slsL1[["DP1.10086.001"]]$"sls_soilMoisture",
                              all_of(joining_cols), moistureSampleID, smSamplingProtocolVersion=samplingProtocolVersion, soilMoisture, smDataQF), by=joining_cols) %>%
      # finally merge with soilpH data
      full_join(dplyr::select(slsL1[["DP1.10086.001"]]$"sls_soilpH",
                              all_of(joining_cols), pHSampleID, pHSamplingProtocolVersion=samplingProtocolVersion, soilInWaterpH, soilInCaClpH, pHDataQF), by=joining_cols)
  } else {
    dat_soil_phys <- NULL
  }

  if(any(grepl('10078', dpID)) & !all(is.na(slsL1[["DP1.10078.001"]]))) {
    dat_soil_chem <- dplyr::select(slsL1[["DP1.10078.001"]]$"sls_soilChemistry",
                                   all_of(joining_cols), cnSampleID, nitrogenPercent, organicCPercent, CNratio, cnTestMethod=testMethod, cnInstrument=instrument, cnDataQF=dataQF)
  } else {
    dat_soil_chem <- NULL
  }

  # TODO: Confirm that the selected columns are sufficient for downstream analysis
  # TODO: Filter by nTransBoutType to remove incubated samples?

  if(!is.null(dat_soil_phys) & !is.null(dat_soil_chem)) {
    message("Returning soil physical and chemical variables at the specified sites and dates.")
    dat_soil <- full_join(dat_soil_phys, dat_soil_chem, by=joining_cols)
  } else if(!is.null(dat_soil_phys) & is.null(dat_soil_chem)) {
    message("Returning soil physical variables (but not chemical variables) at the specified sites and dates.")
    dat_soil <- dat_soil_phys
  } else if(is.null(dat_soil_phys) & !is.null(dat_soil_chem)) {
    message("Returning soil chemical variables (but not physical variables) at the specified sites and dates.")
    dat_soil <- dat_soil_chem
  } else {
    warning("No soil data available at the specified sites and dates. Returning NULL.")
    dat_soil <- NULL
  }

  # download local copy if user provided output dir path
  if(outDir != "") {
    if(!dir.exists(outDir)) {
      dir.create(outDir)
    }
    write.csv(dat_soil, paste0(outDir, "/sls_soilData_", Sys.Date(), ".csv"),
              row.names=F)
    message(paste0("soil data downloaded to: ", outDir, "/sls_soilData_", Sys.Date(), ".csv") )
  }

  return(dat_soil)
}


#' Download Sequence Metadata
#'
#' Loads soil marker gene sequencing metadata for specified target gene, site(s) and date(s),
#' with an option to download output by providing a valid output directory. This function uses
#' \code{\link[neonUtilities]{loadByProduct}} to conduct the downloads.
#'
#' Function by Lee F. Stanish and Clara Qin (2020). Currently available for testing only.
#'
#' @param sites Either the string 'all', meaning all available sites, or a character vector of 4-letter NEON site codes, e.g. c('ONAQ','RMNP'). Defaults to PRESET_SITES parameter in params.R.
#' @param startYrMo,endYrMo Either NA (default), meaning all available dates, or a character vector in the form YYYY-MM, e.g. 2017-01.
#' @param targetGene '16S' or 'ITS'.
#' @param sequencingRuns Either the string 'all', meaning all available sequencing runs, or a character vector of NEON sequencing run IDs, e.g. c('C25G9', 'B69PP').
#' @param dpID NEON data product of interest. Default is soil marker gene sequences, and currently code only works for marker genes data products.
#' @param outDir (Optional) If a local copy of the filtered metadata is desired, provide path to output directory.
#'
#' @return Data frame containing joined records from across the NEON soil marker gene sequence metadata, subsetted according to function arguments.
#' @examples
#' \dontrun{
#' meta <- downloadSequenceMetadataRev('all', '2015-01', '2016-01', '16S') # metadata is not saved to local directory
#' meta <- downloadSequenceMetadataRev('all', '2015-01', '2016-01', '16S', dir='./data/') # metadata is saved to local directory
#' }
downloadSequenceMetadata <- function(sites='all', startYrMo=NA, endYrMo=NA, targetGene= "all",
                                     sequencingRuns = "", dpID = "DP1.10108.001", outDir=PRESET_OUTDIR_SEQMETA) {
  # author: Lee Stanish
  # date: 2020-08-13
  # function loads soil marker gene sequencing metadata for target gene, site(s) and date(s)
  # option to download output by providing a valid output directory
  # sites: character vector of valid site ID's, or 'all' for all sites
  # targetGene: '16S',  'ITS', 'all'
  # startYrMo: start date, format YYYY-MM
  # endYrMo: end date, format YYYY-MM
  # dpID: NEON data product of interest. Default is soil marker gene sequences, and currently code only works for this dpID
  # outDir: directory for outputs. Defaults to output directory in parameters file

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

  # validate output directory
  if(exists("PRESET_OUTDIR_SEQMETA")) {
    outDir <- paste(PRESET_OUTDIR, PRESET_OUTDIR_SEQMETA, sep="/")
  }
  if(!dir.exists(outDir) ) {
    message("Output directory does not exist")
    return(NULL)
  }


  message("loading metadata...")
  mmgL1 <- loadByProduct(dpID, sites, package = 'expanded', check.size = F, startdate = startYrMo, enddate = endYrMo) # output is a list of each metadata file


  # for target data product and targetGene: extract lists into data.frames
  if(grepl("10108", dpID)) {
    seq16S <- mmgL1$mmg_soilMarkerGeneSequencing_16S
    seq16S$targetGene <-"16S rRNA"
    seqITS <- mmgL1$mmg_soilMarkerGeneSequencing_ITS
    seqITS$targetGene <- "ITS"
    pcr16S <- mmgL1$mmg_soilPcrAmplification_16S
    pcrITS <- mmgL1$mmg_soilPcrAmplification_ITS
    raw <- mmgL1$mmg_soilRawDataFiles
    dna <- mmgL1$mmg_soilDnaExtraction
    seq <- rbind(seq16S, seqITS)
    pcr <- rbind(pcr16S, pcrITS)
    varfile <- mmgL1$variables_10108

    if(targetGene=="16S") {
      message("filtering to 16S data")
      seq <- seq16S
      pcr <- pcr16S
    }
    if(targetGene=="ITS") {
      message("filtering to ITS data")
      seq <- seqITS
      pcr <- pcrITS
    }
  }

  if(grepl("20280", dpID)) {
    seq16S <- mmgL1$mmg_benthicMarkerGeneSequencing_16S
    seq16S$targetGene <-"16S rRNA"
    seqITS <- mmgL1$mmg_benthicMarkerGeneSequencing_ITS
    seqITS$targetGene <- "ITS"
    pcr16S <- mmgL1$mmg_benthicPcrAmplification_16S
    pcrITS <- mmgL1$mmg_benthicPcrAmplification_ITS
    raw <- mmgL1$mmg_benthicRawDataFiles
    dna <- mmgL1$mmg_benthicDnaExtraction
    seq <- rbind(seq16S, seqITS)
    pcr <- rbind(pcr16S, pcrITS)
    varfile <- mmgL1$variables_20280

    if(targetGene=="16S") {
      message("filtering to 16S data")
      seq <- seq16S
      pcr <- pcr16S
    }
    if(targetGene=="ITS") {
      message("filtering to ITS data")
      seq <- seqITS
      pcr <- pcrITS
    }
  }

  if(grepl("20282", dpID)) {
    seq16S <- mmgL1$mmg_swMarkerGeneSequencing_16S
    seq16S$targetGene <-"16S rRNA"
    seqITS <- mmgL1$mmg_swMarkerGeneSequencing_ITS
    seqITS$targetGene <- "ITS"
    pcr16S <- mmgL1$mmg_swPcrAmplification_16S
    pcrITS <- mmgL1$mmg_swPcrAmplification_ITS
    raw <- mmgL1$mmg_swRawDataFiles
    dna <- mmgL1$mmg_swDnaExtraction
    seq <- rbind(seq16S, seqITS)
    pcr <- rbind(pcr16S, pcrITS)
    varfile <- mmgL1$variables_20282

    if(targetGene=="16S") {
      message("filtering to 16S data")
      seq <- seq16S
      pcr <- pcr16S
    }
    if(targetGene=="ITS") {
      message("filtering to ITS data")
      seq <- seqITS
      pcr <- pcrITS
    }
  }

  # remove unnecessary/redundant columns from tables
  raw <- dplyr::select(raw, -domainID, -siteID, -namedLocation, -laboratoryName, -sequencingFacilityID, -collectDate, -dnaSampleCode)
  dna <- dplyr::select(dna, -domainID, -siteID, -namedLocation, -laboratoryName, -collectDate)
  pcr <- dplyr::select(pcr, -domainID, -siteID, -namedLocation, -laboratoryName, -collectDate)

  # convert factors to characters (bug in output of loadByProduct)
  i <- sapply(seq, is.factor)
  seq[i] <- lapply(seq[i], as.character)
  j <- sapply(raw, is.factor)
  raw[j] <- lapply(raw[j], as.character)
  j <- sapply(dna, is.factor)
  dna[j] <- lapply(dna[j], as.character)


  # If specified, filter by sequencing run ID
  if(sequencingRuns[1] != "") {
    raw <- raw[which(raw$sequencerRunID %in% sequencingRuns), ]
    # Validate sequencing run ID argument
    if(nrow(raw) == 0) {
      warning("After filtering by specified sequencing run ID(s), no records remain. Double-check your sequencing run  ID(s).")
      return(NULL)
    }
  }

  # Join sequencing metadata with raw data files metadata
  if(targetGene=="16S") {
    if(any(grepl("ITS", raw$rawDataFileName))) {
      rawCleaned <- raw[-grep("ITS", raw$rawDataFileName), ]
    } else {
      rawCleaned <- raw
    }
    joinedTarget <- left_join(rawCleaned, seq, by=c('dnaSampleID', 'sequencerRunID'))
    out <- joinedTarget[!is.na(joinedTarget$uid.y), ]
  }
  if(targetGene=="ITS") {
    if(any(grepl("16S", raw$rawDataFileName))) {
      rawCleaned <- raw[-grep("16S", raw$rawDataFileName), ]
    } else {
      rawCleaned <- raw
    }
    joinedTarget <- left_join(rawCleaned, seq, by=c('dnaSampleID', 'sequencerRunID'))
    out <- joinedTarget[!is.na(joinedTarget$uid.y), ]
  }
  if(targetGene=="all") {
    joinedTarget <- left_join(raw, seq, by=c('dnaSampleID', 'sequencerRunID'))
    out <- joinedTarget[!is.na(joinedTarget$uid.y), ]
    message(paste0(length(grep("16S", out$rawDataFileName)), " 16S records and ", length(grep("ITS", out$rawDataFileName)), " ITS records found."))
  }

  # clean up redundant column names
  names(out) <- gsub("\\.x", ".rawFiles", names(out))
  names(out) <- gsub("\\.y", ".seq", names(out))

  # join with DNA extraction metadata
  outDNA <- left_join(out, dna, by=c('plotID', 'dnaSampleID'))
  # clean up redundant column names
  names(outDNA) <- gsub("\\.x", ".seq", names(outDNA))
  names(outDNA) <- gsub("\\.y", ".dna", names(outDNA))
  names(outDNA)[names(outDNA)=="uid"] <- 'uid.dna'
  names(outDNA)[names(outDNA)=="remarks"] <- 'remarks.dna'
  names(outDNA)[names(outDNA)=="dataQF"] <- 'dataQF.dna'
  names(outDNA)[names(outDNA)=="processedBy"] <- 'processedBy.seq'
  names(outDNA)[names(outDNA)=="processedDate"] <- 'processedDate.dna'
  names(outDNA)[names(outDNA)=="publicationDate"] <- 'publicationDate.dna'
  names(outDNA)[names(outDNA)=="dnaProcessedBy"] <- 'processedBy.dna'

  # join with PCR amplification metadata
  outPCR <- left_join(outDNA, pcr, by=c('plotID', 'dnaSampleID', 'targetGene'))
  names(outPCR)[names(outPCR)=="uid"] <- "uid.pcr"
  names(outPCR)[names(outPCR)=="processedDate"] <- "processedDate.pcr"
  names(outPCR)[names(outPCR)=="testProtocolVersion"] <- "testProtocolVersion.pcr"
  names(outPCR)[names(outPCR)=="qaqcStatus"] <- "qaqcStatus.pcr"
  names(outPCR)[names(outPCR)=="processedBy"] <- "processedBy.pcr"
  names(outPCR)[names(outPCR)=="remarks"] <- "remarks.pcr"
  names(outPCR)[names(outPCR)=="dataQF"] <- "dataQF.pcr"
  names(outPCR)[names(outPCR)=="publicationDate"] <- "publicationDate.pcr"
  names(outPCR)[names(outPCR)=="internalLabID.y"] <- "internalLabID.pcr"

  # download local copy to output dir path
  if(targetGene != "all") {
    write.csv(outPCR, paste0(outDir, "/mmg_soilMetadata_", targetGene, "_", gsub(" |:", "", Sys.time()), ".csv"),
              row.names=FALSE)
  } else {
    out16S <- outPCR[grep("16S", outPCR$targetGene), ]
    outITS <- outPCR[grep("ITS", outPCR$targetGene), ]
    write.csv(out16S, paste0(outDir, "/mmg_soilMetadata_16S_", gsub(" |:", "", Sys.time()), ".csv"),
              row.names=FALSE)
    write.csv(outITS, paste0(outDir, "/mmg_soilMetadata_ITS_", gsub(" |:", "", Sys.time()), ".csv"),
              row.names=FALSE)
  }
  message(paste0("metadata downloaded to: ", outDir) )

  # download variables file (required for zipsByUri)
  write.csv(varfile, paste0(outDir, "/mmg_variables.csv") )
  message(paste0("variables file downloaded to: ", outDir) )

  return(outPCR)
}


#' Trim Primers from 16S Sequences
#'
#' Trims primers from 16S sequences using \code{\link[dada2]{filterAndTrim}}. This function
#' assumes that each read begins with its full primer sequence and operates
#' by truncating the beginning of each read by the length of its primer
#' sequence.
#'
#' @param fn Names of input fastq files, excluding directory path which is specified by dir_in. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param dir_in Directory containing input fastq files.
#' @param dir_out Output directory. If it does not exist, it will be created.
#' @param primer_16S_fwd,primer_16S_rev DNA sequences of 16S forward and reverse primer, respectively
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#' @param post_samplename_pattern1,post_samplename_pattern2 (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#'
#' @return Integer matrix denoting the number of reads remaining after primer-trimming for each input file.
#'
#' @examples
#' \dontrun{
#' trimPrimers16S(c("sample1_R1.fastq", "sample1_R2.fastq", "sample2_R1.fastq", "sample2_R2.fastq"), "path/to/input", "path/to/output", "CCTACGGGNBGCASCAG", "GACTACNVGGGTATCTAATCC", multithread = TRUE)
#' }
trimPrimers16S <- function(fn, dir_in, dir_out, primer_16S_fwd, primer_16S_rev, multithread = MULTITHREAD, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq"#, quiet=T
){
  fn_fullname <- file.path(dir_in, fn)

  fnFs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern1, fn_fullname)]
  fnRs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern2, fn_fullname)]
  if(length(fnFs) + length(fnRs) == 0) warning(paste0("trimPrimer16S: ", "No files found at specified location(s) within ", dir_in, ". Check file path, or post_samplename_pattern argument(s)."))

  # Extract sample names
  sample.names <- getSampleName(fnFs)

  # Create output directory
  if(!dir.exists(dir_out)) dir.create(dir_out)

  # Create primer-trimmed file paths in output directory
  fnFs.cut <- file.path(dir_out, basename(fnFs))
  fnRs.cut <- file.path(dir_out, basename(fnRs))

  # Just trim the length of fwd and rev primers
  qa.out <- filterAndTrim(fnFs, fnFs.cut, fnRs, fnRs.cut, multithread = multithread,
                          matchIDs = TRUE, trimLeft = c(nchar(primer_16S_fwd), nchar(primer_16S_rev)), compress=F)

  colnames(qa.out) <- c("input", "trimmed")
  rownames(qa.out) <- sample.names
  return(qa.out)
}


#' Trim Primers from 16S Sequences (with metadata)
#'
#' Trims primers from 16S sequences using \code{\link[dada2]{filterAndTrim}}. This function
#' assumes that each read begins with its full primer sequence and operates
#' by truncating the beginning of each read by the length of its primer
#' sequence.
#'
#' @param fn Full names of input fastq files, including directory. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param fn_out Full names of output locations, where filtered fastq files will be written.
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param primer_16S_fwd,primer_16S_rev DNA sequences of 16S forward and reverse primer, respectively
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#'
#' @return Integer matrix denoting the number of reads remaining after primer-trimming for each input file.
#'
#' @examples
#' \dontrun{
#' trimPrimers16S(c("sample1_R1.fastq", "sample1_R2.fastq", "sample2_R1.fastq", "sample2_R2.fastq"), "path/to/input", "path/to/output", "CCTACGGGNBGCASCAG", "GACTACNVGGGTATCTAATCC", multithread = TRUE)
#' }
trimPrimers16S2 <- function(fn, fn_out, meta, primer_16S_fwd, primer_16S_rev, multithread = MULTITHREAD) {
  if(length(fn) != length(fn_out)) stop("fn and fn_out must be the same length")

  keep_fn <- file.exists(fn) & !duplicated(fn)
  fn <- fn[keep_fn]
  if(length(fn) == 0) warning(paste0("trimPrimers16S: No files found at specified location(s). Check file paths, or input metadata."))

  # Get metadata matching files
  meta_ext <- matchFastqToMetadata(fn, meta)

  # Reference metadata to retrieve R1 and R2 files
  fn_pairs <- getPairedFastqFiles(fn, meta, value=FALSE)
  fnFs <- fn[fn_pairs[[1]]]
  fnRs <- fn[fn_pairs[[2]]]

  # Confirm target gene
  if(any(!grepl("16S", meta_ext$targetGene))) warning("You are using trimPrimers16S() on some non-16S files. Did you mean to use trimPrimersITS()?")

  # Attempt to get DNA sample IDs to use as row names
  dnaSampleIDs <- meta_ext$dnaSampleID
  fn_as_rownames <- any(is.na(dnaSampleIDs))
  if(fn_as_rownames) {
    warning(sum(is.na(dnaSampleIDs)),
            " files do not have associated dnaSampleIDs in the input metadata. ",
            "Therefore, input R1 file names, instead of dnaSampleIDs, will be used as the rownames for the read tracking table.")
  }

  # Confirm output filenames
  trimFs <- fn_out[keep_fn][fn_pairs[[1]]]
  trimRs <- fn_out[keep_fn][fn_pairs[[2]]]

  # Just trim the length of fwd and rev primers
  qa.out <- filterAndTrim(fnFs, trimFs, fnRs, trimRs, multithread = multithread,
                          matchIDs = TRUE, trimLeft = c(nchar(primer_16S_fwd), nchar(primer_16S_rev)), compress=F)

  colnames(qa.out) <- c("input", "trimmed")
  rownames(qa.out) <- sample.names
  return(qa.out)
}



#' Trim Primers from ITS Sequences
#'
#' Trims primers from ITS sequences using cutadapt. Cutadapt must be installed in order for this to work. Currently only supports R1 (forward-read) files.
#'
#' @param fn Names of input fastq files, excluding directory path which is specified by dir_in. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning. It is assumed that these are R1 (forward-read) files only.
#' @param dir_in Directory containing input fastq files.
#' @param dir_out Output directory. If it does not exist, it will be created.
#' @param primer_ITS_fwd,primer_ITS_rev DNA sequence of the ITS forward and reverse primer, respectively.
#' @param cutadapt_path Default CUTADAPT_PATH in params.R. Path to cutadapt on your file system.
#' @param post_samplename_pattern1,post_samplename_pattern2 (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#' @param very_verbose Default FALSE. Whether to print output from cutadapt. Unlike some other "verbose" arguments associated with the functions in this package, this does not default to VERBOSE in params.R.
#' @param discard_untrimmed Default FALSE. Whether to discard reads where a primer could not be found, leaving only those reads in which a primer has been trimmed.
#'
#' @return No value is returned.
#'
#' @examples
#' \dontrun{
#' trimPrimersITS(c("sample1_ITS_R1.fastq", "sample2_ITS_R1.fastq"), "path/to/input", "path/to/output", "CTTGGTCATTTAGAGGAAGTAA")
#' }
trimPrimersITS <- function(fn, dir_in, dir_out, primer_ITS_fwd, primer_ITS_rev, cutadapt_path = CUTADAPT_PATH, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq", very_verbose=FALSE, discard_untrimmed=FALSE#, quiet=T
){

  fn_fullname <- file.path(dir_in, fn)

  fnFs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern1, fn_fullname)]
  if(length(fnFs) == 0) warning(paste0("trimPrimerITS: ", "No files found at specified location(s) within ", dir_in, ". Check file path, or post_samplename_pattern argument(s)."))

  # Extract sample names
  sample.names <- getSampleName(fnFs)

  # Create output directory
  if(!dir.exists(dir_out)) dir.create(dir_out)

  # Create primer-trimmed file paths in output directory
  fnFs.cut <- file.path(dir_out, basename(fnFs))

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
                                    "--minimum-length", "1", # min length of cutadapted reads: >0
                                    discard_untrimmed_flag,
                                    "-e", "0.2"), # -e 0.2 allows up to 4 mismatched bases
            stdout = ifelse(very_verbose, "", FALSE))
  }
}


#' Trim Primers from ITS Sequences (with metadata)
#'
#' Trims primers from ITS sequences using cutadapt. Cutadapt must be installed in order for this to work. Currently only supports R1 (forward-read) files.
#'
#' @param fn Full names of input fastq files, including directory. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning. It is assumed that these are R1 (forward-read) files only.
#' @param fn_out Full names of output locations, where filtered fastq files will be written.
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param primer_ITS_fwd,primer_ITS_rev DNA sequence of the ITS forward and reverse primer, respectively.
#' @param cutadapt_path Default CUTADAPT_PATH in params.R. Path to cutadapt on your file system.
#' @param very_verbose Default FALSE. Whether to print output from cutadapt. Unlike some other "verbose" arguments associated with the functions in this package, this does not default to VERBOSE in params.R.
#' @param discard_untrimmed Default FALSE. Whether to discard reads where a primer could not be found, leaving only those reads in which a primer has been trimmed.
#'
#' @return No value is returned.
#'
#' @examples
#' \dontrun{
#' trimPrimersITS(c("sample1_ITS_R1.fastq", "sample2_ITS_R1.fastq"), "path/to/input", "path/to/output", "CTTGGTCATTTAGAGGAAGTAA")
#' }
trimPrimersITS2 <- function(fn, fn_out, meta, primer_ITS_fwd, primer_ITS_rev, cutadapt_path = CUTADAPT_PATH, very_verbose=FALSE, discard_untrimmed=FALSE) {
  if(length(fn) != length(fn_out)) stop("fn and fn_out must be the same length")

  keep_fn <- file.exists(fn) & !duplicated(fn)
  fn <- fn[keep_fn]
  if(length(fn) == 0) warning(paste0("trimPrimersITS: No files found at specified location(s). Check file paths, or input metadata."))

  # Reference metadata to retrieve R1 files
  meta_ext <- matchFastqToMetadata(fn, meta)
  fnFs <- fn[grep("R1", meta_ext$rawDataFileDescription)]

  # Confirm target gene
  if(any(!grepl("ITS", meta_ext$targetGene))) warning("You are using trimPrimersITS() on some non-ITS files. Did you mean to use trimPrimers16S()?")

  # Attempt to get DNA sample IDs to use as row names
  dnaSampleIDs <- meta_ext$dnaSampleID
  fn_as_rownames <- any(is.na(dnaSampleIDs))
  if(fn_as_rownames) {
    warning(sum(is.na(dnaSampleIDs)),
            " files do not have associated dnaSampleIDs in the input metadata. ",
            "Therefore, input R1 file names, instead of dnaSampleID, will be used as the rownames for the read tracking table.")
  }

  # Confirm output filenames
  trimFs <- fn_out[keep_fn][grep("R1", meta_ext$rawDataFileDescription)]

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
                                    "-o", trimFs[i], # output files
                                    fnFs[i], # input files
                                    "--minimum-length", "1", # min length of cutadapted reads: >0
                                    discard_untrimmed_flag,
                                    "-e", "0.2"), # -e 0.2 allows up to 4 mismatched bases
            stdout = ifelse(very_verbose, "", FALSE))
  }
}


#' Remove Unmatched Fastq Files (DEPRECATED)
#'
#' NOTE: This function is now deprecated in favor of \code{\link{getPairedFastqFiles}} and \code{\link{removeUnpairedFastqFiles}}.
#'
#' Removes any forward-read files that do not have reverse-read counterparts, and vise versa.
#' This function is necessary because dada2::filterAndTrim() will throw an error
#' if the forward-read files and the reverse-read files are mismatched.
#'
#' @param fnFs Full name(s) of fastq file(s) containing forward-read sequence data.
#' @param fnRs Full name(s) of fastq file(s) containing reverse-read sequence data.
#' @param post_samplename_pattern (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#' @param verbose Whether to print messages regarding which files are matched and which are unmatched.
#'
#' @return List of length 2. The first element is a vector of forward-read files that have reverse-read counterparts; the second element is a vector of reverse-read files that have forward-read counterparts.
#'
#' @examples
#' matched_fn <- remove_unmatched_files(c("sample1_R1.fastq", "sample2_R1.fastq"), c("sample1_R2.fastq", "sample2_R2.fastq", "sample3_R2.fastq"))
#' fnFs <- matched_fn[[1]]
#' fnRs <- matched_fn[[2]]
remove_unmatched_files <- function(fnFs, fnRs, post_samplename_pattern = "_R(1|2).*\\.fastq", verbose=FALSE){
  message("remove_unmatched_files() is deprecated. Please use getPairedFastqFiles() or removeUnpairedFastqFiles() instead.")
  basefilenames_Fs <- sapply(strsplit(fnFs, post_samplename_pattern), `[`, 1)
  basefilenames_Rs <- sapply(strsplit(fnRs, post_samplename_pattern), `[`, 1)
  rm_from_fnFs <- which(!(basefilenames_Fs %in% basefilenames_Rs))
  rm_from_fnRs <- which(!(basefilenames_Rs %in% basefilenames_Fs))

  if (length(c(rm_from_fnFs, rm_from_fnRs)) == 0){
    if(verbose == TRUE) message("All R1/R2 files had counterparts.\n")
  } else {
    for(i in rm_from_fnFs) {
      if(verbose) message(paste(basefilenames_Fs[i], "does not have an R2 counterpart. Omitting from this analysis.\n"))
      fnFs <- fnFs[-rm_from_fnFs]
    }
    for(i in rm_from_fnRs) {
      if(verbose) message(paste(basefilenames_Rs[i], "does not have an R1 counterpart. Omitting from this analysis.\n"))
      fnRs <- fnRs[-rm_from_fnRs]
    }
  }
  return(list(R1=fnFs, R2=fnRs))
}


#' Get Paired Fastq Files
#'
#' Given a vector of fastq file names and corresponding metadata, outputs
#' a list containing matching R1 and R2 files in the same order.
#' Useful for providing files to functions that utilize
#' \code{\link[dada2]{filterAndTrim}}, such as \code{\link{trimPrimer16S}}
#' and \code{\link{runDadaITS}}.
#'
#' See related function: \code{\link{removeUnpairedFastqFiles}}.
#'
#' @param fn Full name(s) of fastq file(s). Tolerant to filenames with appended run ID prefix (e.g. "runB69RN_...") and agnostic to .gz compression.
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param value (Optional) Default TRUE; returns file names. If FALSE, returns indices.
#' @param verbose (Optional) Default TRUE. Whether to print warning messages when metadata does not contain records corresponding to all provided fastq files.
#'
#' @return List of length 2; first element contains matched R1 (forward) file names, and second element contains corresponding R2 (reverse) file names.
#'
#' @examples
#' \dontrun{
#' #' matched_fn <- getPairedFastqFiles(c("sample1_R1.fastq", "sample3_R1.fastq", "sample1_R2.fastq", "sample2_R2.fastq", "sample3_R2.fastq"), meta)
#' fnFs <- matched_fn[[1]] # "sample1_R1.fastq" "sample3_R1.fastq"
#' fnRs <- matched_fn[[2]] # "sample1_R2.fastq" "sample3_R2.fastq"
#' }
getPairedFastqFiles <- function(fn, meta, value=TRUE, verbose=TRUE) {
  meta_ext <- matchFastqToMetadata(fn, meta, verbose=verbose)
  meta_ext$orientation <- if_else(grepl("R1", meta_ext$rawDataFileDescription), "R1",
                                  if_else(grepl("R2", meta_ext$rawDataFileDescription), "R2",
                                          NA_character_))

  meta_ext %>%
    dplyr::filter(!is.na(uid.rawFiles)) %>%
    dplyr::group_by(sequencerRunID, dnaSampleID) %>%
    dplyr::summarise(n_orientations = n_distinct(orientation, na.rm=TRUE), .groups="drop") ->
    meta_summ
  if(!any(meta_summ$n_orientations == 2)) {
    message("There were no matches. Check your inputs and try again.")
    return(NULL)
  } else {
    matched_ids <- meta_summ$dnaSampleID[meta_summ$n_orientations==2]
  }

  if(value==TRUE) {
    return(list(
      meta_ext$file[meta_ext$orientation=="R1" & meta_ext$dnaSampleID %in% matched_ids],
      meta_ext$file[meta_ext$orientation=="R2" & meta_ext$dnaSampleID %in% matched_ids]
    ))
  } else {
    return(list(
      which(meta_ext$orientation=="R1" & meta_ext$dnaSampleID %in% matched_ids),
      which(meta_ext$orientation=="R2" & meta_ext$dnaSampleID %in% matched_ids)
    ))
  }
}


#' Remove Unpaired Fastq Files
#'
#' Given a vector of fastq file names and corresponding metadata, outputs
#' a list containing matching R1 and R2 files in the same order.
#' Useful for providing files to functions that utilize
#' \code{\link[dada2]{filterAndTrim}}, such as \code{\link{trimPrimer16S}}
#' and \code{\link{runDadaITS}}.
#'
#' See related function: \code{\link{getPairedFastqFiles}}.
#'
#' @param fnFs,fnRs Full name(s) of fastq file(s) corresponding to R1 (forward) reads and R2 (reverse) reads, respectively. Tolerant to filenames with appended run ID prefix (e.g. "runB69RN_...") and agnostic to .gz compression.
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param value (Optional) Default TRUE; returns file names. If FALSE, returns indices.
#' @param verbose (Optional) Default TRUE. Whether to print warning messages when metadata does not contain records corresponding to all provided fastq files.
#'
#' @return List of length 2; first element contains matched R1 (forward) file names, and second element contains corresponding R2 (reverse) file names.
#'
#' @examples
#' \dontrun{
#' #' matched_fn <- removeUnpairedFastqFiles(c("sample1_R1.fastq", "sample3_R1.fastq"), c("sample1_R2.fastq", "sample2_R2.fastq", "sample3_R2.fastq"), meta)
#' fnFs <- matched_fn[[1]] # "sample1_R1.fastq" "sample3_R1.fastq"
#' fnRs <- matched_fn[[2]] # "sample1_R2.fastq" "sample3_R2.fastq"
#' }
removeUnpairedFastqFiles <- function(fnFs, fnRs, meta, value=TRUE, verbose=TRUE) {
  meta_ext <- matchFastqToMetadata(c(fnFs, fnRs), meta, verbose=verbose)
  meta_ext$orientation <- c(rep("R1", length(fnFs)), rep("R2", length(fnRs)))

  meta_ext %>%
    dplyr::filter(!is.na(uid.rawFiles)) %>%
    dplyr::group_by(sequencerRunID, dnaSampleID) %>%
    dplyr::summarise(n_orientations = n_distinct(orientation, na.rm=TRUE), .groups="drop") -> meta_summ
  if(!any(meta_summ$n_orientations == 2)) {
    message("There were no matches. Check your inputs and try again.")
    return(NULL)
  } else {
    matched_ids <- meta_summ$dnaSampleID[meta_summ$n_orientations==2]
  }

  if(value==TRUE) {
    return(list(
      meta_ext$file[meta_ext$orientation=="R1" & meta_ext$dnaSampleID %in% matched_ids],
      meta_ext$file[meta_ext$orientation=="R2" & meta_ext$dnaSampleID %in% matched_ids]
    ))
  } else {
    return(list(
      which(meta_ext$orientation=="R1" & meta_ext$dnaSampleID %in% matched_ids),
      which(meta_ext$orientation=="R2" & meta_ext$dnaSampleID %in% matched_ids)
    ))
  }
}


#' Match Fastq Files to Metadata
#'
#' Helper function to match fastq file names to sequence metadata via
#' the sequence metadata table's "rawDataFileName" column. Tolerant to
#' filenames with appended run ID prefix (e.g. "runB69RN_...") and
#' agnostic to .gz compression.
#'
#' Used in \code{\link{getPairedFastqFiles}} and \code{\link{removeUnpairedFastqFiles}}.
#'
#' @param fn Full name(s) of fastq file(s).
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param verbose (Optional) Default TRUE. Whether to print warning messages when metadata does not contain records corresponding to all provided fastq files.
#'
#' @return Data frame, in which the first column contains the fastq file names in the order received, and the remaining columns contain corresponding metadata records.
#'
#' @examples
#' \dontrun{
#' meta_ext <- matchFastqToMetadata(c("sample1_R1.fastq", "sample1_R2.fastq"), )
#' }
matchFastqToMetadata <- function(fn, meta, verbose=TRUE) {
  # Remove runID if appended to beginning of filename
  key <- sub("^run[A-Za-z0-9]*_", "", basename(fn))

  # Append ".gz" to end of filename if missing
  key[!grepl(".gz$", key)] <- paste0(key[!grepl(".gz$", key)], ".gz")

  key_match <- match(key, as.character(meta$rawDataFileName))
  if(any(is.na(key_match))) {
    if(verbose) {
      message("Matching file names to metadata: ", sum(is.na(key_match)), " files did not have matching records in the provided metadata. ",
              "Double-check to ensure that the provided metadata is of the appropriate scope.")
    }
  }
  return(cbind(file = fn, meta[key_match,], stringsAsFactors=FALSE))
}


#' Filter 16S Sequences
#'
#' Applies a quality filter to 16S sequence fastq files via the \code{\link[dada2]{filterAndTrim}} function.
#' It is assumed that both forward- and reverse-read files are included.
#'
#' @param fn Names of input fastq files, excluding directory path which is specified by dir_in. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param dir_in Directory containing input fastq files.
#' @param dir_out Path to output directory where filtered fastq files will be written.
#' @param trunc_qscore Default 23. Quality score at which point to truncate each read, if truncLen is NULL.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#' @param post_samplename_pattern1,post_samplename_pattern2 (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#' @param ... Other arguments to be passed to \code{\link[dada2]{filterAndTrim}}, such as maxEE. See documentation for more details.
#'
#' @return Two-column matrix displaying the number of reads in input vs. output for each file.
qualityFilter16S <- function(fn, dir_in, dir_out, trunc_qscore = 23, multithread = MULTITHREAD, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq", ...){
  fn_fullname <- file.path(dir_in, fn)

  fnFs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern1, fn_fullname)]
  fnRs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern2, fn_fullname)]
  if(length(fnFs) + length(fnRs) == 0) warning(paste0("qualityFilter16S: ", "No files found at specified location(s) within ", dir_in, ". Check file path, or post_samplename_pattern argument(s)."))

  sample.names <- getSampleName(fnFs) # extract sample names
  filtFs <- file.path(dir_out, basename(fnFs)) # create filtered filenames
  filtRs <- file.path(dir_out, basename(fnRs)) # create filtered filenames

  dots <- list(...)

  # If truncLen is NOT in user-provided arguments
  if (!("truncLen" %in% names(dots))) {
    # GETTING TRUNCATION LENGTH USING A SUBSET OF QUALITY SCORES
    n_files <- min(length(unique(sample.names)), 30)

    fwd.trunc.lengths <- list()
    rev.trunc.lengths <- list()
    fwd.trunc.lengths[1:n_files] <- getTruncationLength(fnFs[1:n_files], verbose = VERBOSE,
                                                        qscore = trunc_qscore)
    rev.trunc.lengths[1:n_files] <- getTruncationLength(fnRs[1:n_files], verbose = VERBOSE,
                                                        qscore = trunc_qscore)

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
    if (rev.trunc.length < 200) rev.trunc.length <- 200
    if (fwd.trunc.length < 245) fwd.trunc.length <- 245
    cat(paste0("Fwd truncation length: ", fwd.trunc.length, "\nRev truncation length: ", rev.trunc.length, "\n"))
    truncLen <- c(fwd.trunc.length, rev.trunc.length)

  # Else, if truncLen is in user-provided arguments
  } else {
    message("Using truncLen ", paste(dots$truncLen, collapse=" "))
    arguments <- c(list(fnFs, filtFs, fnRs, filtRs, compress=TRUE, multithread=multithread), dots)
  }

  out <- do.call(filterAndTrim, arguments)

  rownames(out) <- sample.names

  return(out)
}
# Removed this: #' @param truncLen_manual Default NULL. Single integer: truncation length to use across all files. Two-integer vector: truncation length to use for the forward-read and reverse-read files, respectively. If NULL (default), determines truncation length(s) based on \code{\link{getTruncationLength}} with a quality score threshold of trunc_qscore.


#' Filter 16S Sequences (with metadata)
#'
#' Applies a quality filter to 16S sequence fastq files via the \code{\link[dada2]{filterAndTrim}} function.
#' It is assumed that both forward- and reverse-read files are included.
#'
#' @param fn Full names of input fastq files, including directory. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning. It is assumed that these are R1 (forward-read) files only.
#' @param fn_out Full names of output locations, where filtered fastq files will be written.
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param trunc_qscore Default 23. Quality score at which point to truncate each read, if truncLen is NULL.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#' @param ... Other arguments to be passed to \code{\link[dada2]{filterAndTrim}}, such as maxEE and truncLen. See documentation for more details.
#'
#' @return Two-column matrix displaying the number of reads in input vs. output for each file.
qualityFilter16S2 <- function(fn, fn_out, meta, trunc_qscore = 23, multithread = MULTITHREAD, ...){
  if(length(fn) != length(fn_out)) stop("fn and fn_out must be the same length")

  keep_fn <- file.exists(fn) & !duplicated(fn)
  fn <- fn[keep_fn]
  if(length(fn) == 0) warning(paste0("qualityFilter16S: No files found at specified location(s). Check file paths, or input metadata."))

  # Get metadata matching files
  meta_ext <- matchFastqFiles(fn, meta)

  # Reference metadata to retrieve R1 and R2 files
  fn_pairs <- getPairedFastqFiles(fn, meta, value=FALSE)
  fnFs <- fn[fn_pairs[[1]]]
  fnRs <- fn[fn_pairs[[2]]]

  # Confirm target gene
  if(any(!grepl("16S", meta_ext$targetGene))) warning("You are using qualityFilter16S() on some non-16S files. Did you mean to use qualityFilterITS()?")

  # Attempt to get DNA sample IDs to use as row names
  dnaSampleIDs <- meta_ext$dnaSampleID
  fn_as_rownames <- any(is.na(dnaSampleIDs))
  if(fn_as_rownames) {
    warning(sum(is.na(dnaSampleIDs)),
            " files do not have associated dnaSampleIDs in the input metadata. ",
            "Therefore, input R1 file names, instead of dnaSampleIDs, will be used as the rownames for the read tracking table.")
  }

  # Confirm output filenames
  filtFs <- fn_out[keep_fn][fn_pairs[[1]]]
  filtRs <- fn_out[keep_fn][fn_pairs[[2]]]

  # Read in dots args
  dots <- list(...)

  # If truncLen is NOT in user-provided arguments
  if (!("truncLen" %in% names(dots))) {
    # GETTING TRUNCATION LENGTH USING A SUBSET OF QUALITY SCORES
    n_files <- min(length(fnFs), 30)

    fwd.trunc.lengths <- list()
    rev.trunc.lengths <- list()
    fwd.trunc.lengths[1:n_files] <- getTruncationLength(fnFs[1:n_files], verbose = VERBOSE,
                                                        qscore = trunc_qscore)
    rev.trunc.lengths[1:n_files] <- getTruncationLength(fnRs[1:n_files], verbose = VERBOSE,
                                                        qscore = trunc_qscore)

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
    if (rev.trunc.length < 200) rev.trunc.length <- 200
    if (fwd.trunc.length < 245) fwd.trunc.length <- 245
    message("Using auto-selected truncLen:")
    cat(paste0("Fwd truncation length: ", fwd.trunc.length, "\nRev truncation length: ", rev.trunc.length, "\n"))
    dots$truncLen <- c(fwd.trunc.length, rev.trunc.length)

    # Else, if truncLen is in user-provided arguments
  } else {
    message("Using user-provided truncLen: ", paste(dots$truncLen, collapse=" "))
  }

  # Run dada::filterAndTrim
  arguments <- c(list(fnFs, filtFs, fnRs, filtRs, compress=TRUE, multithread=multithread), dots)
  out <- do.call(filterAndTrim, arguments)

  if(fn_as_rownames) {
    rownames(out) <- fnFs
  } else {
    rownames(out) <- dnaSampleIDs
  }

  return(out)
}


#' Filter ITS Sequences
#'
#' Applies a quality filter to ITS sequence fastq files via the \code{\link[dada2]{filterAndTrim}} function.
#' Currently only supports filtering forward-read sequences.
#'
#' @param fn Names of input fastq files, excluding directory path which is specified by dir_in. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning. It is assumed that these are R1 (forward-read) files only.
#' @param dir_in Directory containing input fastq files.
#' @param dir_out Path to output directory where filtered fastq files will be written.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#' @param post_samplename_pattern1,post_samplename_pattern2 (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#' @param ... Other arguments to be passed to \code{\link[dada2]{filterAndTrim}}, such as maxEE. See documentation for more details.
#'
#' @return Two-column matrix displaying the number of reads in input vs. output for each file.
qualityFilterITS <- function(fn, dir_in, dir_out, multithread = MULTITHREAD, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq", ...){

  fn_fullname <- file.path(dir_in, fn)

  fnFs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern1, fn_fullname)]
  if(length(fnFs) == 0) warning(paste0("trimPrimerITS: ", "No files found at specified location(s) within ", dir_in, ". Check file path, or post_samplename_pattern argument(s)."))

  sample.names <- getSampleName(fnFs) # extract sample names

  gz_ind <- grepl("\\.gz$", basename(fnFs))
  if(!all(gz_ind)) {
    filtFs <- file.path(dir_out, basename(fnFs))
    filtFs[!gz_ind] <- paste0(filtFs, ".gz")
  } else {
    filtFs <- file.path(dir_out, basename(fnFs))
  }

  out <- filterAndTrim(fnFs, filtFs,
                       compress = TRUE,
                       matchIDs = TRUE,
                       multithread = multithread,
                       ...)

  rownames(out) <- sample.names

  return(out)
}

#' Filter ITS Sequences (with metadata)
#'
#' Applies a quality filter to ITS sequence fastq files via the \code{\link[dada2]{filterAndTrim}} function.
#' Currently only supports filtering forward-read sequences.
#'
#' @param fn Full names of input fastq files, including directory. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning. It is assumed that these are R1 (forward-read) files only.
#' @param fn_out Full names of output locations, where filtered fastq files will be written.
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#' @param ... Other arguments to be passed to \code{\link[dada2]{filterAndTrim}}, such as maxEE. See documentation for more details.
#'
#' @return Two-column matrix displaying the number of reads in input vs. output for each file.
qualityFilterITS2 <- function(fn, fn_out, meta, multithread = MULTITHREAD, ...){
  if(length(fn) != length(fn_out)) stop("fn and fn_out must be the same length")

  keep_fn <- file.exists(fn) & !duplicated(fn)
  fn <- fn[keep_fn]
  if(length(fn) == 0) warning(paste0("qualityFilterITS: No files found at specified location(s). Check file paths, or input metadata."))

  # Reference metadata to retrieve R1 files
  meta_ext <- matchFastqToMetadata(fn, meta)
  fnFs <- fn[grep("R1", meta_ext$rawDataFileDescription)]

  # Confirm target gene
  if(any(!grepl("ITS", meta_ext$targetGene))) warning("You are using qualityFilterITS() on some non-ITS files. Did you mean to use qualityFilter16S()?")

  # Attempt to get DNA sample IDs to use as row names
  dnaSampleIDs <- meta_ext$dnaSampleID
  fn_as_rownames <- any(is.na(dnaSampleIDs))
  if(fn_as_rownames) {
    warning(sum(is.na(dnaSampleIDs)),
            " files do not have associated dnaSampleIDs in the input metadata. ",
            "Therefore, input R1 file names, instead of dnaSampleID, will be used as the rownames for the read tracking table.")
  }

  # Confirm output filenames
  filtFs <- fn_out[keep_fn][grep("R1", meta_ext$rawDataFileDescription)]

  # Run quality filter
  out <- filterAndTrim(fnFs, filtFs,
                       compress = TRUE,
                       multithread = multithread,
                       ...)

  if(fn_as_rownames) {
    rownames(out) <- fnFs
  } else {
    rownames(out) <- dnaSampleIDs
  }
  return(out)
}



#' Get Truncation Length
#'
#' Decides on truncation length for trimmed 16S reads based on quality score means. Default cutoff is a score of 30. Warns you if the beginning (first 10) bases are low-quality, but returns the first low-quality base after the 10th base. If no bases are below score, returns the last base. The truncation lengths can be aggregated (e.g. minimum) and used as an argument in \code{\link[dada2]{filterAndTrim}}.
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
#' @param fn Names of input fastq files, excluding directory path which is specified by dir_in. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param dir_in Directory containing input fastq files.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading.
#' @param VERBOSE Default FALSE. Whether to print messages regarding the dimensions of the resulting sequence table and the distribution of sequence lengths.
#' @param seed (Optional) Integer to use as random seed for reproducibility.
#' @param nbases (Optional) Number of bases to use for learning errors. Default 1e7.
#'
#' @return A list of three elements. \strong{seqtab} is the sequence table before removing chimeras, \strong{seqtab.nochim} is the sequence table after removing chimeras, and \strong{track} is a data frame displaying the number of reads remaining for each sample at various points throughout the processing pipeline.
#'
#' @examples
#' \dontrun{
#' seqtab.list <- runDada16S(c("sample1_R1.fastq", "sample1_R2.fastq", "sample2_R1.fastq", "sample2_R2.fastq"), seed=1010100)
#' }
runDada16S <- function(fn, dir_in, multithread = MULTITHREAD, verbose = FALSE, seed = NULL, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq", nbases=1e7){
  if (!is.null(seed)) set.seed(seed)

  fn_fullname <- file.path(dir_in, fn)

  filtFs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern1, fn_fullname)]
  filtRs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern2, fn_fullname)]
  if(length(filtFs) + length(filtRs) == 0) warning(paste0("runDada16S: ", "No files found at specified location(s) within ", dir_in, ". Check file path, or post_samplename_pattern argument(s)."))

  # Keep files with counterpart
  matched_files <- remove_unmatched_files(filtFs, filtRs)
  filtFs <- matched_files[["R1"]]
  filtRs <- matched_files[["R2"]]

  # Create outnames
  sample.names <- sapply(strsplit(basename(filtFs), paste0("(",post_samplename_pattern1,")|(",post_samplename_pattern2,")")), `[`, 1)
  sample.namesR <- sapply(strsplit(basename(filtRs), paste0("(",post_samplename_pattern1,")|(",post_samplename_pattern2,")")), `[`, 1)
  # sample.names <- sapply(strsplit(basename(filtFs), "_R(1|2).fastq"), `[`, 1) # Assumes filename = samplename_RX.fastq.gz
  # sample.namesR <- sapply(strsplit(basename(filtRs), "_R(1|2).fastq"), `[`, 1) # Assumes filename = samplename_RX.fastq.gz
  if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  # Learn forward and reverse error rates
  errF <- learnErrors(filtFs, nbases=nbases, randomize=TRUE, multithread=multithread, verbose=verbose)
  errR <- learnErrors(filtRs, nbases=nbases, randomize=TRUE, multithread=multithread, verbose=verbose)

  # Create output vectors
  derepF <- derepR <- ddF <- ddR <- mergers <- vector("list", length(sample.names))
  names(derepF) <- names(derepR) <- names(ddF) <- names(ddR) <- names(mergers) <- sample.names

  # Sample inference and merger of paired-end reads
  for(i in 1:length(sample.names)) {
    sam <- sample.names[[i]]
    cat("Processing:", sam, "\n")

    derepF[[i]] <- derepFastq(filtFs[[i]])
    derepR[[i]] <- derepFastq(filtRs[[i]])

    ddF[[i]] <- dada(derepF[[i]], err=errF, multithread=multithread, verbose=verbose)
    ddR[[i]] <- dada(derepR[[i]], err=errR, multithread=multithread, verbose=verbose)
    mergers[[i]] <- mergePairs(ddF[[i]], derepF[[i]], ddR[[i]], derepR[[i]], maxMismatch=1, minOverlap = 6, verbose=verbose)
    cat(paste("Exact sequence variants inferred for sample:", sam,". \n"))
  }

  # Construct sequence table and remove chimeras
  seqtab <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=multithread, verbose=verbose)

  if(verbose){
    cat(paste0("\n\nDimensions of ESV table: ", dim(seqtab)[1], " samples, ", dim(seqtab)[2], " ESVs\n"))
    cat("\nDistribution of sequence lengths (should be bimodal for v3v4 region):")
    print(table(nchar(getSequences(seqtab))))
    cat(paste0("\n", round(sum(seqtab.nochim)/sum(seqtab), 3)*100, "% reads remain after removal of chimeras"))
  }

  # Track reads through this part of the pipeline
  suppressWarnings({
    track <- Reduce(
    function(x, y, ...) transform(merge(x, y, by=0, all = TRUE, ...), row.names=Row.names, Row.names=NULL),
    list(sapply(derepF, getN),
         sapply(derepR, getN),
         sapply(ddF, getN),
         denoisedR = sapply(ddR, getN),
         sapply(mergers, getN),
         rowSums(seqtab.nochim))
    )
  })
  names(track) <- c("derepF", "derepR", "denoisedF", "denoisedR", "merged", "nonchim")
  track[is.na(track)] <- 0

  return(list("seqtab" = seqtab,
              "seqtab.nochim" = seqtab.nochim,
              "track" = track))
}


#' Run Dada on Paired-End 16S Sequences (with metadata)
#'
#' Runs the Dada algorithm to infer sample composition from paired-end 16S fastq files.
#' This implementation is based on Ben Callahan's vignette at \url{https://benjjneb.github.io/dada2/bigdata_paired.html}.
#'
#' @param fn Full names of input fastq files, including directory. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param out_seqtab,out_track (Optional) File locations where sequence table and read-tracking table (respectively) will be saved as csv files. If blank (default), will not save to file.
#' @param remove_chimeras (Optional) Default TRUE. Whether to remove chimeras from the sequence table. Currently only supports removal using \code{\link[dada2]{removeBimeraDenovo}} with the consensus method.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading.
#' @param verbose Default FALSE. Whether to print messages regarding the samples being processed, dimensions of the resulting sequence table, and the distribution of sequence lengths.
#' @param seed (Optional) Integer to use as random seed for reproducibility.
#' @param nbases (Optional) Number of bases to use for learning errors. Default 1e7.
#'
#' @return A list of two elements. \strong{seqtab} is the sequence table, with chimeras removed if remove_chimeras == TRUE, and \strong{track} is a data frame displaying the number of reads remaining for each sample at various points throughout the processing pipeline.
#'
#' @examples
#' \dontrun{
#' seqtab.list <- runDada16S(c("sample1_R1.fastq", "sample1_R2.fastq", "sample2_R1.fastq", "sample2_R2.fastq"), seed=1010100)
#' }
runDada16S2 <- function(fn, meta, out_seqtab = "", out_track = "", remove_chimeras = TRUE, multithread = MULTITHREAD, verbose = FALSE, seed = NULL, nbases=1e7){
  if (!is.null(seed)) set.seed(seed)

  keep_fn <- file.exists(fn) & !duplicated(fn)
  fn <- fn[keep_fn]
  if(length(fn) == 0) warning(paste0("runDada16S: No files found at specified location(s). Check file paths, or input metadata."))

  # Get metadata matching files
  meta_ext <- matchFastqToMetadata(fn, meta)

  # Reference metadata to retrieve R1 and R2 files
  fn_pairs <- getPairedFastqFiles(fn, meta, value=FALSE)
  fnFs <- fn[fn_pairs[[1]]]
  fnRs <- fn[fn_pairs[[2]]]

  # Confirm target gene
  if(any(!grepl("16S", meta_ext$targetGene))) warning("You are using runDada16S() on some non-16S files. Did you mean to use runDadaITS()?")

  # Attempt to get DNA sample IDs to use as row names
  dnaSampleIDs <- as.character(meta_ext$dnaSampleID[match(fnFs, meta_ext$file)])
  fn_as_rownames <- any(is.na(dnaSampleIDs))
  if(fn_as_rownames) {
    warning(sum(is.na(dnaSampleIDs)),
            " files do not have associated dnaSampleIDs in the input metadata. ",
            "Therefore, input R1 file names, instead of dnaSampleIDs, will be used as the rownames for the read tracking table.")
  }

  names(fnFs) <- dnaSampleIDs
  names(fnRs) <- dnaSampleIDs

  # Learn forward and reverse error rates
  errF <- learnErrors(fnFs, nbases=nbases, randomize=TRUE, multithread=multithread, verbose=verbose)
  errR <- learnErrors(fnRs, nbases=nbases, randomize=TRUE, multithread=multithread, verbose=verbose)

  # Create output vectors
  derepF <- derepR <- ddF <- ddR <- mergers <- vector("list", length(dnaSampleIDs))
  names(derepF) <- names(derepR) <- names(ddF) <- names(ddR) <- names(mergers) <- dnaSampleIDs

  # Sample inference and merger of paired-end reads
  message("runDada16S progress out of ", length(dnaSampleIDs), " samples:")
  progressbar <- txtProgressBar(min = 0, max = length(dnaSampleIDs), style = 3)
  for(i in 1:length(dnaSampleIDs)) {
    sam <- dnaSampleIDs[[i]]
    if(verbose) cat("Processing:", sam, "\n")

    derepF[[i]] <- derepFastq(fnFs[[i]])
    derepR[[i]] <- derepFastq(fnRs[[i]])

    ddF[[i]] <- dada(derepF[[i]], err=errF, multithread=multithread, verbose=verbose)
    ddR[[i]] <- dada(derepR[[i]], err=errR, multithread=multithread, verbose=verbose)
    if(verbose) {
      mergers[[i]] <- mergePairs(ddF[[i]], derepF[[i]], ddR[[i]], derepR[[i]], maxMismatch=1, minOverlap = 6, verbose=TRUE)
    } else {
      mergers[[i]] <- suppressMessages(mergePairs(ddF[[i]], derepF[[i]], ddR[[i]], derepR[[i]], maxMismatch=1, minOverlap = 6, verbose=FALSE))
    }
    if(verbose) cat(paste("Exact sequence variants inferred for sample:", sam,". \n"))
    setTxtProgressBar(progressbar, i)
  }
  close(progressbar)

  # Construct sequence table and remove chimeras
  if(verbose) {
    seqtab <- makeSequenceTable(mergers)
  } else {
    seqtab <- suppressMessages(makeSequenceTable(mergers))
  }
  if(remove_chimeras) seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=multithread, verbose=verbose)

  if(verbose){
    cat(paste0("\nDimensions of ESV table: ", dim(seqtab)[1], " samples, ", dim(seqtab)[2], " ESVs\n"))
    cat("\nDistribution of sequence lengths (should be bimodal for v3v4 region):")
    print(table(nchar(getSequences(seqtab))))
    cat(paste0("\n", round(sum(seqtab.nochim)/sum(seqtab), 3)*100, "% reads remain after removal of chimeras"))
  }

  # Track reads through this part of the pipeline
  suppressWarnings({
    track <- Reduce(
      function(x, y, ...) transform(merge(x, y, by=0, all = TRUE, ...), row.names=Row.names, Row.names=NULL),
      list(sapply(derepF, getN),
           sapply(derepR, getN),
           sapply(ddF, getN),
           denoisedR = sapply(ddR, getN),
           suppressMessages(sapply(mergers, getN)),
           rowSums(seqtab.nochim))
    )
  })
  names(track) <- c("derepF", "derepR", "denoisedF", "denoisedR", "merged", "nonchim")
  track[is.na(track)] <- 0

  if(out_seqtab != "") {
    if(!dir.exists(dirname(out_seqtab))) {
      dir.create(dirname(out_seqtab), recursive=TRUE)
      message("Output directory created for sequence table: ", dirname(out_seqtab))
    }
    if(remove_chimeras) {
      write.csv(seqtab.nochim, out_seqtab)
    } else {
      write.csv(seqtab, out_seqtab)
    }
  }
  if(out_track != "") {
    if(!dir.exists(dirname(out_track))) {
      dir.create(dirname(out_track), recursive=TRUE)
      message("Output directory created for read tracking table: ", dirname(out_track))
    }
    write.csv(track, out_track)
  }

  if(remove_chimeras) {
    return(list("seqtab" = seqtab.nochim,
                "track" = track))
  } else {
    return(list("seqtab" = seqtab,
                "track" = track))
  }
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
#' @param fn Names of input fastq files, excluding directory path which is specified by dir_in. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning. It is assumed that these are R1 (forward-read) files only.
#' @param dir_in Directory containing input fastq files.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading.
#' @param verbose Default FALSE. Whether to print messages regarding the dereplication step, the denoising step, and the dimensions of the resulting sequence table and the distribution of sequence lengths.
#' @param seed (Optional) Integer to use as random seed for reproducibility.
#' @param post_samplename_pattern1,post_samplename_pattern2 (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#' @param nbases (Optional) Number of bases to use for learning errors. Default 1e7.
#'
#' @return A list of three elements. \strong{seqtab} is the sequence table before removing chimeras, \strong{seqtab.nochim} is the sequence table after removing chimeras, and \strong{track} is a data frame displaying the number of reads remaining for each sample at various points throughout the processing pipeline.
#'
#' @examples
#' \dontrun{
#' seqtab.list <- runDadaITS(c("sample1_R1.fastq", "sample1_R2.fastq", "sample2_R1.fastq", "sample2_R2.fastq"), './seq/filtered/', seed=1010100)
#' }
runDadaITS <- function(fn, dir_in, multithread = MULTITHREAD, verbose = FALSE, seed = NULL, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq", nbases = 1e7){
  if (!is.null(seed)) set.seed(seed)

  fn_fullname <- file.path(dir_in, fn)

  filtFs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern1, fn_fullname)]
  if(length(filtFs) == 0) warning(paste0("runDadaITS: ", "No files found at specified location(s) within ", dir_in, ". Check file path, or post_samplename_pattern argument(s)."))

  # Create outnames
  sample.names <- sapply(strsplit(basename(filtFs), paste0("(",post_samplename_pattern1,")|(",post_samplename_pattern2,")")), `[`, 1)
  names(filtFs) <- sample.names

  # Learn error rates
  errF <- learnErrors(filtFs, nbases=nbases, randomize=TRUE, multithread=multithread, verbose=verbose)

  # Create output vectors
  derepF <- ddF <- vector("list", length(sample.names))
  names(derepF) <- names(ddF) <- sample.names

  # Sample inference and merger of paired-end reads
  for(i in 1:length(sample.names)) {
    sam <- sample.names[[i]]
    cat("Processing:", sam, "\n")
    derepF[[i]] <- derepFastq(filtFs[[i]], verbose=verbose)
    ddF[[i]] <- dada(derepF[[i]], err=errF, multithread=multithread, verbose=verbose)
    cat(paste("Exact sequence variants inferred for sample:", sam,". \n"))
  }
  # rm(derepF);

  # Construct sequence table and remove chimeras
  seqtab <- makeSequenceTable(ddF)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=multithread, verbose=verbose)

  if(verbose){
    cat(paste0("\n\nDimensions of ESV table: ", dim(seqtab)[1], " samples, ", dim(seqtab)[2], " ESVs\n"))
    cat("\nDistribution of sequence lengths:")
    print(table(nchar(getSequences(seqtab))))
    cat(paste0("\n", round(sum(seqtab.nochim)/sum(seqtab), 3)*100, "% reads remain after removal of chimeras"))
  }

  # Track reads through this part of the pipeline
  track <- Reduce(
    function(x, y, ...) transform(merge(x, y, by=0, all = TRUE, ...), row.names=Row.names, Row.names=NULL),
    list(sapply(derepF, getN),
         sapply(ddF, getN),
         rowSums(seqtab.nochim))
  )
  names(track) <- c("derepF", "denoisedF", "nonchim")
  track[is.na(track)] <- 0

  return(list("seqtab" = seqtab,
              "seqtab.nochim" = seqtab.nochim,
              "track" = track))
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


#' Get Sample Name from Fastq Filename (DEPRECATED)
#'
#' Deprecated in favor of code that reads dnaSampleID from sequence metadata.
#'
#' Extracts the sample names from fastq filenames, given regular expressions
#' denoting the different endings of R1 and R2 filenames.
#'
#' @param fn Fastq filename(s).
#' @param post_samplename_pattern1,post_samplename_pattern2 Default "_R1.*\\.fastq" and "_R2.*\\.fastq". Regular expressions denoting the ends of the R1 and R2 filenames, respectively, which will be trimmed away, leaving identical names between corresponding R1 and R2 files.
#'
#' @return Character vector; sample names of the input fastq filenames.
#'
#' @examples
getSampleName <- function(fn, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq") {
  message("getSampleName() is deprecated. Please retrieve dnaSampleID from sequence metadata instead.")
  sapply(strsplit(basename(fn), paste0("(",post_samplename_pattern1,")|(",post_samplename_pattern2, ")")), `[`, 1)
}


#' Parse Parameter Values from Rownames
#'
#' Parses parameter values from rownames. Particularly useful for sensitivity analyses, when rows have been constructed to contain information about the parameter values, e.g. "maxEE.R_8_truncLen.R_250".
#'
#' @param df Data frame with rownames containing parameter values of the form [param_name1]_[param_value1]_[param_name2]_[param_value2].
#' @param keep_rownames (Default FALSE) Whether to keep the rownames of the original dataframe.
#' @param keep_original_cols (Default TRUE) Whether to keep the columns of the original dataframe.
#' @param param1,param2 Parameter names embedded within the rownames of df.
#'
#' @return A data frame containing columns corresponding to the values of the parameters of interest, extracted from the rownames of the original data frame.
#'
#' @examples \dontrun{ trackReadsWithParamVals <- parseParamsFromRownames(trackReads, PARAM1, PARAM2) }
parseParamsFromRownames <- function(df, param1, param2, keep_rownames=FALSE, keep_original_cols=TRUE) {
  if(class(df)=="matrix") df <- as.data.frame(df)
  param1.0 <- as.numeric(unlist(regmatches(rownames(df), gregexpr(paste0(param1, "_\\K[0-9]*"), rownames(df), perl=TRUE))))
  param2.0 <- as.numeric(unlist(regmatches(rownames(df), gregexpr(paste0(param2, "_\\K[0-9]*"), rownames(df), perl=TRUE))))
  runID0 <- unlist(regmatches(rownames(df), gregexpr("\\|\\Krun[A-Za-z0-9]*", rownames(df), perl=TRUE)))
  sampleID0 <- unlist(regmatches(rownames(df), gregexpr("\\|\\K.*", rownames(df), perl=TRUE)))
  if(length(param1.0)==nrow(df)) df[[param1]] <- param1.0
  if(length(param2.0)==nrow(df)) df[[param2]] <- param2.0
  if(length(runID0)==nrow(df)) df[["runID"]] <- runID0
  if(length(sampleID0)==nrow(df)) df[["sampleID"]] <- sampleID0
  if(!keep_rownames) {
    rownames(df) <- NULL
  }
  if(keep_original_cols) {
    return(df)
  } else {
    param1 <- enquo(param1)
    param2 <- enquo(param2)
    return(dplyr::select(df, any_of(c(!!param1, !!param2, "runID", "sampleID"))))
  }
}



#' Plot Expected Errors Profile
#'
#' Plots the cumulative expected errors of a fastq file or collection of fastq files over sequence positions. This function extends \code{\link[dada2]{plotQualityProfile}} and makes use of its code.
#'
#' @param fn Full name(s) of input fastq file(s), including directory.
#' @param n (Optional) Default 5e+5. The number of reads to sample when processing fastq files; passed to \code{\link[ShortRead]{qa}}.
#' @param aggregate (Optiona) Default FALSE. If inputting multiple fastq files, whether to aggregate their expected errors into one summary plot. If "aggregate" and "include_quality_profile" are both TRUE, also aggregates the quality profiles into one summary plot.
#' @param logEE (Optional) Default TRUE. Whether to log10-transform the Y-axis of the expected error profile.
#' @param include_quality_profile (Optional) Default FALSE. Whether to include the quality profile in the return value (see below).
#'
#' @return If include_quality_profile = FALSE (default), returns a ggplot of the expected errors. If include_quality_profile = TRUE, returns a length of list two: (1) a ggplot of the expected errors, and (2) a ggplot of the quality profile(s).
#'
#' @examples
#' \dontrun{
#' fnFs <- c("sample1_R1.fastq", "sample2_R1.fastq", "sample3_R1.fastq")
#' plotEEProfile(fnFs) # produces a multi-faceted ggplot
#' plotEEProfile(fnFs, aggregate=TRUE) # produces a single-faceted ggplot
#' p <- plotEEProfile(fnFs, aggregate=TRUE, include_quality_profile=TRUE) # produces a list containing two single-faceted ggplots:
#' plot(p[[1]]) # expected errors profile
#' plot(p[[2]]) # quality profile
#' }
plotEEProfile <- function(fn, n=5e+5, aggregate=FALSE, logEE=TRUE, include_quality_profile=FALSE) {
  if(length(fn) > 20 & missing(aggregate)) {
    message("More files than recommended for individual plotting. Setting aggregate = TRUE. ",
            "Disable by explicitly setting aggregate = FALSE.")
    aggregate <- TRUE
  }
  statdf <- data.frame(Cycle = integer(0), Mean = numeric(0),
                       Q25 = numeric(0), Q50 = numeric(0), Q75 = numeric(0),
                       Cum = numeric(0), file = character(0))
  anndf <- data.frame(minScore = numeric(0), label = character(0),
                      rclabel = character(0), rc = numeric(0), file = character(0))
  FIRST <- TRUE
  message("Calculating summary statistics for plotting...")
  progressbar <- txtProgressBar(min = 0, max = length(fn[!is.na(fn)]), style = 3)
  for (i in 1:length(fn[!is.na(fn)])) {
    f <- fn[!is.na(fn)][i]
    srqa <- qa(f, n = n)
    df <- srqa[["perCycle"]]$quality
    rc <- sum(srqa[["readCounts"]]$read)
    if (rc >= n) {
      rclabel <- paste("Reads >= ", n)
    } else {
      rclabel <- paste("Reads: ", rc)
    }
    means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count,
                                                          df$Cycle)
    get_quant <- function(xx, yy, q) {
      xx[which(cumsum(yy)/sum(yy) >= q)][[1]]
    }
    q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score,
                                                     foo$Count, 0.25), simplify = TRUE)
    q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score,
                                                     foo$Count, 0.5), simplify = TRUE)
    q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score,
                                                     foo$Count, 0.75), simplify = TRUE)
    cums <- by(df, df$Cycle, function(foo) sum(foo$Count),
               simplify = TRUE)
    if (!all(sapply(list(names(q25s), names(q50s), names(q75s),
                         names(cums)), identical, rownames(means)))) {
      stop("Calculated quantiles/means weren't compatible.")
    }
    if (FIRST) {
      plotdf <- cbind(df, file = basename(f))
      FIRST <- FALSE
    } else {
      plotdf <- rbind(plotdf, cbind(df, file = basename(f)))
    }
    n.min <- min(n, rc)
    statdf <- rbind(statdf, data.frame(Cycle = as.integer(rownames(means)),
                                       Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s),
                                       Q75 = as.vector(q75s), Cum = 10 * as.vector(cums)/n.min, #instead of /rc
                                       file = basename(f)))
    anndf <- rbind(anndf, data.frame(minScore = min(df$Score),
                                     label = basename(f), rclabel = rclabel, rc = rc,
                                     file = basename(f)))
    setTxtProgressBar(progressbar, i)
  }
  close(progressbar)

  anndf$minScore <- min(anndf$minScore)

  ee_label0 <- data.frame(ee=2^seq(0,5), label=paste0("maxEE = ", 2^seq(0,5)))

  # If aggregating files
  if (aggregate) {
    plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf,
                                sum)
    plotdf.summary$label <- paste(nrow(anndf), "files (aggregated)")
    means <- rowsum(plotdf.summary$Score * plotdf.summary$Count,
                    plotdf.summary$Cycle)/rowsum(plotdf.summary$Count,
                                                 plotdf.summary$Cycle)
    q25s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score,
                                                                             foo$Count, 0.25), simplify = TRUE)
    q50s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score,
                                                                             foo$Count, 0.5), simplify = TRUE)
    q75s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score,
                                                                             foo$Count, 0.75), simplify = TRUE)
    cums <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) sum(foo$Count),
               simplify = TRUE)
    n.min <- ifelse(anndf$rc < n, anndf$rc, n)
    n.agg <- mean(n.min)
    statdf.summary <- data.frame(Cycle = as.integer(rownames(means)),
                                 Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s),
                                 Q75 = as.vector(q75s), Cum = 10 * as.vector(cums)/n.agg)

    if(include_quality_profile==TRUE) {
      p_q <- ggplot(data = plotdf.summary, aes(x = Cycle, y = Score)) +
      geom_tile(aes(fill = Count)) + scale_fill_gradient(low = "#F5F5F5", high = "black") +
      geom_line(data = statdf.summary, aes(y = Mean), color = "#66C2A5") +
      geom_line(data = statdf.summary, aes(y = Q25), color = "#FC8D62", size = 0.25, linetype = "dashed") +
      geom_line(data = statdf.summary, aes(y = Q50), color = "#FC8D62", size = 0.25) +
      geom_line(data = statdf.summary, aes(y = Q75), color = "#FC8D62", size = 0.25, linetype = "dashed") +
      ylab("Quality Score") +
      xlab("Cycle") +
      annotate("text", x = 0, y = 0, label = sprintf("Total reads: %d", sum(anndf$rc)), color = "red", hjust = 0) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      guides(fill = FALSE) +
      facet_wrap(~label) +
      ylim(c(0, NA))
      if (length(unique(statdf$Cum)) > 1) {
        p_q <- p_q +
          geom_line(data = statdf.summary, aes(y = Cum/nrow(anndf)),
                    color = "red", size = 0.25, linetype = "solid") +
          scale_y_continuous(sec.axis = sec_axis(~. * 10, breaks = c(0, 100), labels = c("0%", "100%"))) +
          theme(axis.text.y.right = element_text(color = "red"),
                axis.title.y.right = element_text(color = "red"))
      }
    }

    eedf <- dplyr::mutate(statdf.summary, across(c(starts_with("Q"), Mean), cumulativeExpectedErrors))
    eedf <- tidyr::pivot_longer(eedf, c(starts_with("Q"), Mean), names_to="series", values_to="expected_errors")

    p_ee <- ggplot(eedf) +
      geom_line(aes(x=Cycle, y=expected_errors, col=series, linetype=series, size=series)) +
      scale_color_manual(values=c("#66C2A5", "#FC8D62", "#FC8D62", "#FC8D62")) +
      scale_linetype_manual(values=c("solid", "dashed", "solid", "dashed")) +
      scale_size_manual(values=c(1, 0.5, 0.5, 0.5)) +
      labs(y="Cum. expected errors", col="") +
      guides(col=FALSE, linetype=FALSE, size=FALSE)

  # If not aggregating files
  } else {
    if(include_quality_profile) {
      p_q <- ggplot(data = plotdf, aes(x = Cycle, y = Score)) +
        geom_tile(aes(fill = Count)) +
        scale_fill_gradient(low = "#F5F5F5", high = "black") +
        geom_line(data = statdf, aes(y = Mean), color = "#66C2A5") +
        geom_line(data = statdf, aes(y = Q25), color = "#FC8D62", size = 0.25, linetype = "dashed") +
        geom_line(data = statdf, aes(y = Q50), color = "#FC8D62", size = 0.25) +
        geom_line(data = statdf, aes(y = Q75), color = "#FC8D62", size = 0.25, linetype = "dashed") +
        ylab("Quality Score") + xlab("Cycle") + theme_bw() +
        theme(panel.grid = element_blank()) + guides(fill = FALSE) +
        geom_text(data = anndf, aes(x = 0, label = rclabel,
                                    y = 0), color = "red", hjust = 0) +
        facet_wrap(~file) +
        ylim(c(0, NA))
      if (length(unique(statdf$Cum)) > 1) {
        p_q <- p_q + geom_line(data = statdf, aes(y = Cum), color = "red", size = 0.25, linetype = "solid") +
          scale_y_continuous(sec.axis = sec_axis(~. * 10, breaks = c(0, 100), labels = c("0%", "100%"))) +
          theme(axis.text.y.right = element_text(color = "red"),
                axis.title.y.right = element_text(color = "red"))
      }
    }

    statdf %>%
      tidyr::pivot_longer(c(starts_with("Q"), Mean), names_to="series", values_to="quality") %>%
      group_by(file, series) %>%
      arrange(Cycle) %>%
      dplyr::mutate(expected_errors = cumulativeExpectedErrors(quality)) ->
      eedf

    ee_label <- ee_label0[ee_label0$ee < max(eedf$expected_errors),]

    p_ee <- ggplot(eedf) +
      geom_line(aes(x=Cycle, y=expected_errors, col=series, linetype=series, size=series)) +
      scale_color_manual(values=c("#66C2A5", "#FC8D62", "#FC8D62", "#FC8D62")) +
      scale_linetype_manual(values=c("solid", "dashed", "solid", "dashed")) +
      scale_size_manual(values=c(1, 0.5, 0.5, 0.5)) +
      labs(y="Cum. expected errors", col="") +
      facet_wrap(~file) +
      guides(col=FALSE, linetype=FALSE, size=FALSE)
  }

  ee_label <- ee_label0[ee_label0$ee < max(eedf$expected_errors),]
  if(nrow(ee_label) > 0) {
    p_ee <- p_ee +
      geom_text(data=ee_label, aes(x=0, y=ee, label=label), vjust=-0.2, hjust=0, col="grey50", size=3) +
      geom_hline(data=ee_label, aes(yintercept=ee), lwd=0.5, col="grey50")
  }

  if(logEE) {
    p_ee <- p_ee + scale_y_continuous(trans="log10")
  }

  if(include_quality_profile) {
    return(list(p_ee, p_q))
  } else {
    return(p_ee)
  }
}

cumulativeExpectedErrors <- function(Q) {
  cumsum(10^(-Q/10))
}
