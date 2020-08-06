# Utilities for use in other R scripts in NEON/DoB microbial community analysis

# This version of utils.R began as a copy of ./code/utils.R, and is intended for
# writing documentation of new functions using Roxygen. It does not include all
# of the functions in the original ./code/utils.R, like get.sample.name().
# Those functions are either copied from or closely adapted from the DADA tutorials,
# so they should be included in the vignettes, but not claimed as a new function.

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
#' write_sample_subset_params(writeDir="./raw_data/", sites="all", startYrMo="2016-01", endYrMo="2019-01", target_genes="ITS", sequencing_runs=c("B69PP","C25G9"))
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
#' uses the neonUtilities package to conduct the downloads.
#'
#' @param sites Either the string 'all', meaning all available sites, or a character vector of 4-letter NEON site codes, e.g. c('ONAQ','RMNP'). Defaults to PRESET_SITES parameter in params.R.
#' @param startYrMo,endYrMo Either NA, meaning all available dates, or a character vector in the form YYYY-MM, e.g. 2017-01. Defaults to PRESET_START_YR_MO in params.R.
#' @param outdir Location where output files are saved. Defaults to PRESET_OUTDIR_SEQMETA in params.R.
#' @param checkFileSize TRUE or FALSE. Whether the user should be told the total file size before downloading. Defaults to PRESET_CHECK_FILE_SIZE in params.R.
#' @param return_data Whether to return part of the metadata: the table NEON table "mmg_soilRawDataFiles.csv" within NEON data product DP1.10108. If FALSE, simply downloads data to outdir and returns no value. Defaults to PRESET_RETURN_DATA.
#' @param target_genes 'ITS', '16S', or 'all'. Defaults to TARGET_GENE in params.R.
#' @param sequencing_runs Either the string 'all', meaning all available sequencing runs, or a character vector of NEON sequencing run IDs, e.g. c('C25G9', 'B69PP'). Defaults to SEQUENCING_RUNS parameter in params.R.
#' @param overwrite TRUE or FALSE. If there is a previous download of this metadata in outdir, whether to overwrite that download.
#'
#' @return If return_data==TRUE, returns the mmg_soilRawDataFiles table from NEON.DP1.10108 ("Soil microbe marker gene sequences"). Otherwise, no value is returned.
downloadSequenceMetadata <- function(sites = PRESET_SITES, startYrMo = PRESET_START_YR_MO, endYrMo = PRESET_END_YR_MO,
                                     outdir = PRESET_OUTDIR_SEQMETA, checkFileSize = PRESET_CHECK_FILE_SIZE, return_data = PRESET_RETURN_DATA,
                                     target_genes = TARGET_GENE, sequencing_runs = SEQUENCING_RUNS,
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


#' Download NEON Marker Gene Sequencing Raw Data
#'
#' Downloads NEON raw sequence data files to the specified filepath, by referencing
#' URLs in the NEON table "mmg_soilRawDataFiles.csv" within data product DP1.10108
#' ("Soil microbe marker gene sequences"). If the latter has not been downloaded,
#' calls downloadSequenceMetadata().
#'
#' @param sites Either the string 'all', meaning all available sites, or a character vector of 4-letter NEON site codes, e.g. c('ONAQ','RMNP'). Defaults to PRESET_SITES parameter in params.R.
#' @param startYrMo,endYrMo Either NA, meaning all available dates, or a character vector in the form YYYY-MM, e.g. 2017-01. Defaults to PRESET_START_YR_MO in params.R.
#' @param outdir Location where output files are saved. Defaults to PRESET_OUTDIR_SEQMETA in params.R.
#' @param checkFileSize TRUE or FALSE. Whether the user should be told the total file size before downloading. Defaults to PRESET_CHECK_FILE_SIZE in params.R.
#' @param return_data Whether to return metadata. If FALSE, simply downloads data to outdir and returns no value. Defaults to PRESET_RETURN_DATA.
#' @param target_genes 'ITS', '16S', or 'all'. Defaults to TARGET_GENE in params.R.
#' @param sequencing_runs Either the string 'all', meaning all available sequencing runs, or a character vector of NEON sequencing run IDs, e.g. c('C25G9', 'B69PP'). Defaults to SEQUENCING_RUNS parameter in params.R.
#' @param overwrite TRUE or FALSE. If there is a previous download of sequence metadata in outdir, whether to overwrite that download.
#'
#' @return If return_data==TRUE, returns the mmg_soilRawDataFiles table from NEON.DP1.10108 ("Soil microbe marker gene sequences"). Otherwise, no value is returned.
downloadRawSequenceData <- function(sites = PRESET_SITES, startYrMo = PRESET_START_YR_MO, endYrMo = PRESET_END_YR_MO,
                                    outdir = PRESET_OUTDIR_SEQUENCE, checkFileSize = PRESET_CHECK_FILE_SIZE, return_data = PRESET_RETURN_DATA,
                                    target_genes = TARGET_GENE, sequencing_runs = SEQUENCING_RUNS,
                                    overwrite = FALSE) {

  library(utils)
  library(dplyr)

  metadata <- downloadSequenceMetadata(sites, startYrMo, endYrMo, checkFileSize=FALSE, return_data=TRUE,
                                       overwrite=overwrite)

  # If necessary, subset to download specific sequencing runs
  if(length(sequencing_runs) > 1) {
    if(!any(sequencing_runs=="all")) {
      metadata <- dplyr::filter(metadata, sequencerRunID %in% sequencing_runs)
    }
  }

  u.urls <- unique(metadata$rawDataFilePath)
  fileNms <- gsub('^.*\\/', "", u.urls)

  # If necessary, subset to download only ITS or 16S
  if(length(target_genes)==1) {
    if(target_genes=="ITS") {
      keep_ind <- grep("_ITS_",fileNms)
      fileNms <- fileNms[keep_ind]
      u.urls <- u.urls[keep_ind]
    }
    if(target_genes=="16S") {
      keep_ind <- grep("_16S_",fileNms)
      fileNms <- fileNms[keep_ind]
      u.urls <- u.urls[keep_ind]
    }
  } else {
    stop("Error: This implementation only allows 'TARGET_GENE' parameters of length 1")
  }

  print(paste("There are", length(u.urls), "unique (zipped) raw sequence files to download.") )

  for(i in 1:length(u.urls)) {
    download.file(url=as.character(u.urls[i]), destfile = ifelse(dir.exists(outdir),
                                                                 paste(outdir, fileNms[i], sep="/"),
                                                                 paste(getwd(), fileNms[i], sep="/" )) )
    if(dir.exists(outdir)) {
      print(paste("Finished downloading", paste(outdir, fileNms[i], sep="/")))
    } else {
      print(paste("Finished downloading", paste(getwd(), fileNms[i], sep="/" )))
    }
  }

  if(return_data) {
    return(metadata)
  }
}


#' Download NEON Soil Data Associated with Marker Gene Sequencing Data
#'
#' Downloads the following data products to the specified filepath:
#' (1) DP1.10078: "Soil chemical properties (Distributed periodic)";
#' (2) DP1.10086: "Soil physical properties (Distributed periodic)".
#' This function uses the neonUtilities package to conduct the downloads.
#'
#' @param sites Either the string 'all', meaning all available sites, or a character vector of 4-letter NEON site codes, e.g. c('ONAQ','RMNP'). Defaults to PRESET_SITES parameter in params.R.
#' @param startYrMo,endYrMo Either NA, meaning all available dates, or a character vector in the form YYYY-MM, e.g. 2017-01. Defaults to PRESET_START_YR_MO in params.R.
#' @param outdir Location where output files are saved. Defaults to PRESET_OUTDIR_SEQMETA in params.R.
#' @param checkFileSize TRUE or FALSE. Whether the user should be told the total file size before downloading. Defaults to PRESET_CHECK_FILE_SIZE in params.R.
#' @param return_data Whether to return metadata. If FALSE, simply downloads data to outdir and returns no value. Defaults to PRESET_RETURN_DATA.
#' @param overwrite TRUE or FALSE. If there is a previous download of this metadata in outdir, whether to overwrite that download.
#'
#' @return If return_data==TRUE, returns a dataframe consisting of joined soil data records from DP1.10078 ("Soil chemical properties (Distributed periodic)") and DP1.10086 ("Soil physical properties (Distributed periodic)"). Otherwise, no value is returned.
downloadRawSoilData <- function(sites = PRESET_SITES, startYrMo = PRESET_START_YR_MO, endYrMo = PRESET_END_YR_MO,
                                outdir = PRESET_OUTDIR_SOIL, checkFileSize = PRESET_CHECK_FILE_SIZE, return_data = PRESET_RETURN_DATA,
                                overwrite = FALSE) {
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


# NEXT FUNCTION CREATED BY LEE F STANISH. TESTING USE ONLY FOR NOW #
#### FUNCTION downloadSequenceMetadataRev ####
downloadSequenceMetadataRev <- function(sites='all', startYrMo, endYrMo, targetGene= "",
                                        sequencingRuns = "", dpID = "DP1.10108.001", dir="") {
  # author: Lee Stanish
  # date: 2020-06-16
  # function loads soil marker gene sequencing metadata for target gene, site(s) and date(s)
  # option to download output by providing a valid output directory
  # sites: character vector of valid site ID's, or 'all' for all sites
  # targetGene: '16S' or 'ITS'
  # startYrMo: start date, format YYYY-MM
  # endYrMo: end date, format YYYY-MM
  # dpID: NEON data product of interest. Default is soil marker gene sequences, and currently code only works for this dpID
  # dir (optional): If a local copy of the filtered metadata is desired, provide path to output dir

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
  if(!grepl("16S|ITS", targetGene)) {
    message("Invalid targetGene: must be either '16S' or 'ITS'")
    return(NULL)
  } else {
    targetGene <- targetGene
  }

  # How can we validate sequencing runs? Some match the 5-character format, while others are like "150925_M02149_0250_000000000-AHCGJ"
  # # validate sequencing runs
  # if(!any(grepl("[A-Z0-9]{5}|^all$", sequencingRuns))) {
  #   message("Invalid sequencingRuns: must be 'all' or follow convention of 5-character run IDs, e.g. 'C25G9'.")
  #   return(NULL)
  # } else {
  #   sequencingRuns <- sequencingRuns
  # }

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
    if(targetGene=="16S") {
      message("filtering to 16S data")
      seq <- mmgL1$mmg_soilMarkerGeneSequencing_16S
      raw <- mmgL1$mmg_soilRawDataFiles

    } else {
      message("filtering to ITS data")
      seq <- mmgL1$mmg_soilMarkerGeneSequencing_ITS
      raw <- mmgL1$mmg_soilRawDataFiles
    }
  }

  if(grepl("20280", dpID)) {
    if(targetGene=="16S") {
      message("filtering to 16S data")
      seq <- mmgL1$mmg_benthicMarkerGeneSequencing_16S
      raw <- mmgL1$mmg_benthicRawDataFiles
    } else {
      message("filtering to ITS data")
      seq <- mmgL1$mmg_benthicMarkerGeneSequencing_ITS
      raw <- mmgL1$mmg_benthicRawDataFiles
    }
  }

  if(grepl("20282", dpID)) {
    if(targetGene=="16S") {
      message("filtering to 16S data")
      seq <- mmgL1$mmg_swMarkerGeneSequencing_16S
      raw <- mmgL1$mmg_swRawDataFiles
    } else {
      message("filtering to ITS data")
      seq <- mmgL1$mmg_swMarkerGeneSequencing_ITS
      raw <- mmgL1$mmg_swRawDataFiles
    }
  }

  # convert factors to characters (bug in output of loadByProduct)
  i <- sapply(seq, is.factor)
  seq[i] <- lapply(seq[i], as.character)
  j <- sapply(raw, is.factor)
  raw[j] <- lapply(raw[j], as.character)

  # If specified, filter by sequencing run ID
  if(sequencingRuns!="") {
    if(!("all" %in% sequencingRuns)) {
      if(any(!(raw$sequencerRunID %in% sequencingRuns))) {
        raw <- raw[-which(!(raw$sequencerRunID %in% sequencingRuns)), ]
      }
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
  } else if(targetGene=="ITS") {
    if(any(grepl("16S", raw$rawDataFileName))) {
      rawCleaned <- raw[-grep("16S", raw$rawDataFileName), ]
    } else {
      rawCleaned <- raw
    }
    joinedTarget <- left_join(rawCleaned, seq, by=c('dnaSampleID', 'sequencerRunID', 'internalLabID'))
    out <- joinedTarget[!is.na(joinedTarget$uid.y), ]
  }

  # download local copy if user provided output dir path
  if(dir != "") {
    if(!dir.exists(dir)) {
      dir.create(dir)
    }
    write.csv(out, paste0(dir, "/mmg_soilMetadata_", targetGene, "_", Sys.Date(), ".csv"),
              row.names=F)
    print(paste0("metadata downloaded to: ", dir, "/mmg_soilMetadata_", targetGene, "_", Sys.Date(), ".csv") )
  }
  return(out)

  ### END FUNCTION ###
}


meta <- downloadSequenceMetadataRev(sites="TALL", startYrMo="2014-01", endYrMo="2015-01", targetGene="ITS", dir="/data/ZHULAB/NEON_DOB/sequence_metadata")
