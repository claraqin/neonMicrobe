#' Download NEON Marker Gene Sequencing Raw Data
#'
#' Downloads NEON raw sequence data files to the specified filepath, by referencing
#' URLs in the metadata output from \code{\link{downloadSequenceMetadata}}.
#'
#' @param metadata The output of downloadSequenceMetadata(). Must be provided as either the data.frame returned by downloadSequenceMetadata() or as a filepath to the csv file produced by downloadSequenceMetadata() when outdir is provided.
#' @param outdir Location where output files are saved. Defaults to PRESET_OUTDIR_SEQUENCE in params.R.
#' @param ignore_tar_files If TRUE (default), does not download tar files. Each tar file is a batch containing an entire sequence run of fastq files. The tar file structure will soon be deprecated.
#' @param verbose If TRUE, prints status messages and progress bars associated with file downloads.
#'
#' @return Returns (invisibly) a list of integer codes: 0 indicates success of downloads and a non-zero integer indicates failure. See the help page for \code{\link[utils]{download.file}} for more details.
downloadRawSequenceData <- function(metadata, outDir = PRESET_OUTDIR_SEQUENCE,
                                    ignore_tar_files=TRUE, checkSize=TRUE, verbose=FALSE) {

  library(utils)
  options(stringsAsFactors = FALSE)

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

  # Get unique file names
  metadata.u <- metadata[!duplicated(metadata$rawDataFilePath), ]
  print(paste("There are", nrow(metadata.u), "unique raw sequence files to download."))

  #Loop to check existence and cumulative size of files
  cat("checking file sizes...\n")
  fileSize <- 0
  idx <- 1
  idxrem <- vector()
  for(i in 1:nrow(metadata.u)) {
    # get file metadata
    response <- httr::HEAD(metadata.u$rawDataFilePath[i])
    # check for file found
    if(is.null(httr::headers(response)[["Content-Length"]])) {
      cat(paste('No file found for url ', metadata.u$rawDataFilePath[i], '. Skipping\n', sep=''))
      idxrem[idx] <- c(idxrem, i)
      idx <- idx + 1
    } else {
      # grab file size
      fileSize <- fileSize + as.numeric(httr::headers(response)[["Content-Length"]])
    }
  }
  #  Sum up file sizes and convert bytes to MB
  totalFileSize <- fileSize/1e6

  # Remove missing files from download list
  if(length(idxrem)>0) {
    metadata.u <- metadata.u[-idxrem, ]
  }

  if(checkSize==TRUE) {
    resp <- readline(paste("Continuing will download",nrow(metadata.u), "files totaling approximately",
                           totalFileSize, "MB. Do you want to proceed? y/n: ", sep=" "))
    if(!(resp %in% c("y","Y"))) stop("Stopping")
  }else{
    cat("Downloading",length(idx), "files totaling approximately",totalFileSize," MB.\n")
  }


  download_success <- list()
  for(i in 1:nrow(metadata.u)) {
    tryCatch({
      download_success[[i]] <- download.file(
        url = metadata.u$rawDataFilePath[i],
        destfile = ifelse(dir.exists(outDir), paste(outDir, metadata.u$rawDataFileName[i], sep="/"), paste(getwd(), metadata.u$rawDataFileName[i], sep="/" ) ),
        quiet = !verbose)
    }, error = function(e) { # Occasionally an error arises because _fastq should be replaced by .fastq
      tryCatch({
        revised_url <- sub("_fastq", ".fastq", as.character(metadata.u$rawDataFilePath[i]))
        download_success[[i]] <- download.file(
          url = revised_url,
          destfile = ifelse(dir.exists(outDir), paste(outDir, metadata.u$rawDataFileName[i], sep="/"), paste(getwd(), metadata.u$rawDataFileName[i], sep="/" )),
          quiet = !verbose)
      }, error = function(f) {
        message("Could not download from URL: ", metadata.u$rawDataFilePath[i])
        download_success[[i]] <- 2
      })
    })
    if(dir.exists(outDir)) {
      if(verbose) message("Finished downloading ", paste(outDir, metadata.u$rawDataFileName[i], sep="/"), ".\n")
    } else {
      if(verbose) message("Finished downloading ", paste(getwd(), metadata.u$rawDataFileName[i], sep="/" ), ".\n")
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
#' @param metadata The output of downloadSequenceMetadata(). Must be provided as either the data.frame returned by downloadSequenceMetadata() or as a filepath to the csv file produced by downloadSequenceMetadata() when outDir is provided.
#' @param outdir_sequence Default PRESET_OUTDIR_SEQUENCE. Directory where raw sequence files can be found before reorganizing.
#' @param verbose If TRUE, prints message each time a file is reorganized.
#'
#' @return Character vector of the files (including files within tarballs) that were successfully reorganized. If no files were successfully reorganized, returns no value.
organizeRawSequenceData <- function(fn, metadata, outdir_sequence = PRESET_OUTDIR_SEQUENCE, verbose = TRUE) {
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
#' @param sites Either the string 'all' (default), meaning all available sites, or a character vector of 4-letter NEON site codes, e.g. c('ONAQ','RMNP'). Defaults to PRESET_SITES parameter in params.R.
#' @param startYrMo,endYrMo Either NA (default), meaning all available dates, or a character vector in the form YYYY-MM, e.g. 2017-01. Defaults to PRESET_START_YR_MO in params.R.
#' @param dpID NEON data product(s) of interest. Default is both DP1.10078.001 ("Soil chemical properties (Distributed periodic)") and DP1.10086.001 ("Soil physical properties (Distributed periodic)").
#' @param outDir Default PRESET_OUTDIR_SOIL If a local copy of the filtered metadata is desired, provide path to output directory.
#'
#' @return If return_data==TRUE, returns a dataframe consisting of joined soil data records from DP1.10078 ("Soil chemical properties (Distributed periodic)") and DP1.10086 ("Soil physical properties (Distributed periodic)"). Otherwise, no value is returned.
downloadRawSoilData <- function(sites='all', startYrMo = NA, endYrMo = NA,
                                dpID = c("DP1.10078.001", "DP1.10086.001"), outDir=PRESET_OUTDIR_SOIL) {

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
    outDir <- PRESET_OUTDIR_SEQMETA
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

#' QC Sequence Metadata
#'
#' Performs basic QAQC checks on sequence metadata prior to downloading sequence data and performing bioinformatics processing.
#' Running this function will remove metadata records for samples that do not meet user specifications. This will reduce the number of sequence files that are downloaded to only those that will be used for analysis, thereby saving file space and reducing download times.
#'
#' @param metadata The output of downloadSequenceMetadata(). Must be provided as either the data.frame returned by downloadSequenceMetadata() or as a filepath to the csv file produced by downloadSequenceMetadata() when outdir is provided.
#' @param outDir Default PRESET_OUTDIR_SEQMETA. Directory where raw sequence files can be found before reorganizing.
#' @param pairedReads "Y" (default) or "N". Should the forward reads for a sample be removed if the corresponding reverse read is missing? If "Y", then only samples that have both the forward (R1) and reverse (R2) reads will be retained.
#' @param rmDupes TRUE (default) or FALSE. Should records with duplicated dnaSampleIDs be removed? If TRUE, then only the first records encountered for a particular dnaSampleID will be retained.
#'
#' @return QC'd dataframe is returned as an object and saved as csv file.
qcMetadata <- function(metadata, outDir=PRESET_OUTDIR_SEQMETA, pairedReads="Y", rmDupes=TRUE, rmFlagged="N", verbose=FALSE) {
  library(plyr)
  options(stringsAsFactors = FALSE)

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

  # validate pairedReads
  if(!(pairedReads %in% c("Y", "N")) ) {
    stop("value for argument pairedReads invalid. Must be 'Y' or 'N'.")
  }

  # validate rmFlagged
  if(!(rmFlagged %in% c("Y", "N")) ) {
    stop("value for argument rmFlagged invalid. Must be 'Y' or 'N'.")
  }

  # get targetGene and confirm that only one targetGene is in input data set
  targetGene <- unique(metadata$targetGene)
  if(length(targetGene) > 1) {
    stop("more than one targetGene in input data set. Only one targetGene can be QCed at a time.")
  }

  # create output folder for QCed metadata #
  qcDir <- paste0(outDir, '/QC_output/')
  if(!dir.exists(qcDir) ) {
    dir.create(qcDir)
  }

  # Print size of dataset
  print(paste("Input dataset contains", nrow(metadata), "rows.") )


  # Remove flagged records, if rmFlagged="Y"
  if(rmFlagged=="Y") {
    print("Removing records with existing quality flag...")
    # Define flag values for removal based on LOV values #
    flagVals <- "Fail|legacyData"
    # Remove NA values #
    metadata[is.na(metadata)] <- ""
    flagFields <- grep("qaqcStatus|dataQF", names(metadata))
    ind <- vector()
    for(i in flagFields) {
      flagged <- grep(flagVals, metadata[,i])
      ind <- c(ind,flagged)
    }
    if(length(ind)==0) {
      print("No flagged records found.")
    } else {
      ind <- unique(ind)
      numDupes <- length(ind)
      print(paste0(length(ind), " flagged records found. Removing flagged record(s).") )
      metadata <- metadata[-ind, ]
    }
    Sys.sleep(3)
  }

  # check for and remove duplicate sequence file names
  print("QC check for duplicate sequence file names...")
  # Add pause #
  Sys.sleep(1)
  dupeSeqIDs <- as.character(metadata$rawDataFileName[duplicated(metadata$rawDataFileName)] )
  if(length(dupeSeqIDs)==0) {
    print("QC check Pass. No duplicate sequence file names.")
  } else {
    numDupes <- length(dupeSeqIDs)
    print(paste0("QC check Fail. ", numDupes, " duplicate sequence file names found. Removing duplicated file(s).") )
    if(verbose) print(paste0("Removing duplicated row: ", which(duplicated(metadata$rawDataFileName)) ) )
    metadata <- metadata[!duplicated(metadata$rawDataFileName), ]
  }


  # check for and flag duplicate dnaSampleIDs
  print("QC checking duplicate dnaSampleIDs...")
  # Add pause #
  Sys.sleep(2)
  metadata$runDir <- ""
  metadata$runDir[grep("R1", metadata$rawDataFileDescription)] <- "R1"
  metadata$runDir[grep("R2", metadata$rawDataFileDescription)] <- "R2"
  dnaIDsPerRunDir <- paste(metadata$dnaSampleID, metadata$sequencerRunID, metadata$runDir, sep="-")
  metadata$duplicateDnaSampleIDFlag <- "0"
  if(any(duplicated(dnaIDsPerRunDir)) ) {
    ind <- which(duplicated(dnaIDsPerRunDir))
    metadata$duplicateDnaSampleIDFlag[ind] <- "1"
    if(rmDupes==TRUE) {
      metadata <- metadata[-ind,]
      print("Duplicated dnaSampleID record(s) found and removed.")
    } else {
      print("Duplicated dnaSampleID(s) found. Flagging affected record(s).")
    }
  } else {
    print("No duplicated dnaSampleID records found.")
  }

  # Subset qced data to un-flagged records
  metaFlagged <- metadata[metadata$duplicateDnaSampleIDFlag=="1", ]
  metaNotFlagged <- metadata[metadata$duplicateDnaSampleIDFlag=="0", ]

  metaNotFlagged$dnaSampleID <- as.character(metaNotFlagged$dnaSampleID)

  # Check existence of R1 (and R2 based on user input) #
  dnaSampTab <- data.frame(with(metaNotFlagged, table(dnaSampleID, runDir)) )

  # convert factors to characters
  dfType <- sapply(dnaSampTab, class)
  colsToFix <- names(dnaSampTab[which(dfType=='factor')])
  dnaSampTab[colsToFix] <- sapply(dnaSampTab[colsToFix], as.character)

  # Handle sequence data with missing run direction
  missingR1 <- dnaSampTab$dnaSampleID[which(dnaSampTab$Freq[dnaSampTab$runDir=="R1"]==0)]
  missingR2 <- dnaSampTab$dnaSampleID[which(dnaSampTab$Freq[dnaSampTab$runDir=="R2"]==0)]
  metaNotFlagged$runDirFlag <- "0"
  print("Check for missing forward or reverse read...")
  # Add pause #
  Sys.sleep(2)
  # Handle sequence data with missing R2 data
  if(length(missingR2)>0) {
    print(paste0("Reverse read missing from ", length(missingR2), " records") )
    # If specified, remove R1 file
    if(pairedReads=="Y") {
      print("Removing R1 files lacking a matching R2 file (default action when pairedReads='Y')" )
      metaNotFlagged <- metaNotFlagged[-intersect(which(metaNotFlagged$runDir=="R1"), which(metaNotFlagged$dnaSampleID %in% missingR2) ), ]
    } else {
      metaNotFlagged$runDirFlag[intersect(which(metaNotFlagged$runDir=="R1"), which(metaNotFlagged$dnaSampleID %in% missingR2) )] <- "1"
      print("Flagging R1 files lacking a matching R2 file (default action when pairedReads='N')" )
    }
  }
  # Handle sequence data with missing R1 data
  if(length(missingR1)>0) {
    print(paste0("Forward read missing from ", length(missingR1), " record(s). Removing R2 file for affected sequence data set(s).") )
    # remove R2 file
    metaNotFlagged <- metaNotFlagged[-intersect(which(metaNotFlagged$runDir=="R2"), which(metaNotFlagged$dnaSampleID %in% missingR1) ), ]
  }

  # Recombine original flagged records and remaining records post-initial flagging.
  out <- suppressMessages(plyr::join(metaNotFlagged, metaFlagged))

  write.csv(out, paste0(qcDir, "mmg_metadata_", gsub("\\ ", "", targetGene), "_QCed_", gsub("-", "", Sys.Date()), '.csv'), row.names = FALSE)
  cat(paste("Output QCed file contains", nrow(out), "rows. File saved to the following directory:", qcDir))
  cat("\nNOTE: Always review output before proceeding with analysis.")
  return(out)
}
