# Utilities for use in other R scripts in NEON/DoB microbial community analysis

# This line is an addition for example.

# This line is an addition for test.

library(neonUtilities)
library(utils)
library(dplyr)

# Load parameters from params.R
source("./code/params.R")

## Function to write site and date range parameters. To be called within other functions.
#
write_sample_subset_params <- function(writeDir, sites, startYrMo, endYrMo,
                                       target_genes="", sequencing_runs="") {
  # paste("writing in ", writeDir)
  write.table(
    data.frame(x1=c("sites", "startYrMo", "endYrMo", "target_genes", "sequencing_runs"),
               x2=c(paste0(sites, collapse=","), startYrMo, endYrMo, paste0(target_genes, collapse=","), paste0(sequencing_runs, collapse=","))),
    file = paste(writeDir, SAMPLE_SUBSET_PARAMS_FILENAME, sep="/"),
    sep=":", col.names=FALSE, quote=FALSE, row.names=FALSE
  )
}
## END FUNCTION


## Function issues warning that metadata has already been downloaded, and indicates the sites
## and date ranges for which it was downloaded. To be called within other functions.
# TODO: May want to consider an "Overwrite? (y/n)" option instead of requiring
#       user to manually overwrite if they want to redownload.
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
## END FUNCTION


## Function downloads the metadata for NEON marker gene sequencing data products 
#
# Metadata is just one table withing DP1.10108 ("Soil microbe marker gene sequences")
# called "mmg_soilRawDataFiles.csv"
downloadSequenceMetadata <- function(sites = PRESET_SITES, startYrMo = PRESET_START_YR_MO, endYrMo = PRESET_END_YR_MO, 
                                     outdir = PRESET_OUTDIR_SEQMETA, checkFileSize = PRESET_CHECK_FILE_SIZE, return_data = PRESET_RETURN_DATA,
                                     target_genes = TARGET_GENE, sequencing_runs = SEQUENCING_RUNS,
                                     overwrite = FALSE) {
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
## END FUNCTION ##


## Function downloads the metadata for NEON marker gene sequencing data products 
## AND downloads the NEON raw sequence data files
downloadRawSequenceData <- function(sites = PRESET_SITES, startYrMo = PRESET_START_YR_MO, endYrMo = PRESET_END_YR_MO, 
                                    outdir = PRESET_OUTDIR_SEQUENCE, checkFileSize = PRESET_CHECK_FILE_SIZE, return_data = PRESET_RETURN_DATA,
                                    target_genes = TARGET_GENE, sequencing_runs = SEQUENCING_RUNS,
                                    overwrite = FALSE) {
  
  metadata <- downloadSequenceMetadata(sites, startYrMo, endYrMo, checkFileSize=FALSE, return_data=TRUE,
                                       overwrite=overwrite)
  
  
  # If necessary, subset to download specific sequencing runs
  if(length(sequencing_runs) > 1) {
    if(!any(sequencing_runs=="all")) {
      metadata <- filter(metadata, sequencerRunID %in% sequencing_runs)
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
## END FUNCTION ##


## Function downloads NEON soil data to be associated with soil microbial
## community composition data
#
# Data product are
# - DP1.10078: "Soil chemical properties (Distributed periodic)"
# - DP1.10086: "Soil physical properties (Distributed periodic)"
downloadRawSoilData <- function(sites = PRESET_SITES, startYrMo = PRESET_START_YR_MO, endYrMo = PRESET_END_YR_MO, 
                                outdir = PRESET_OUTDIR_SOIL, checkFileSize = PRESET_CHECK_FILE_SIZE, return_data = PRESET_RETURN_DATA,
                                overwrite = FALSE) {
  # TODO: Why does this function sometimes return a few warnings of the form:
  # 1: In UseMethod("depth") :
  # no applicable method for 'depth' applied to an object of class "NULL"
  
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
## END FUNCTION ##


## Function downloads ALL NEON soil data, using preset parameters
#
# Wrapper for downloadRawSoilData,
# parameterized for socs-stats.ucsc.edu server,
# and will pull data up to current month.
downloadAllRawSoilData <- function() {
  return(downloadRawSoilData(PRESET_SITES, PRESET_START_YR_MO, PRESET_END_YR_MO, 
                             PRESET_OUTDIR_SOIL, PRESET_CHECK_FILE_SIZE, PRESET_RETURN_DATA))
}
## END FUNCTION ##


## Function to get all orients of primers for DADA2 workflow
## (From DADA2 ITS Pipeline Workflow 1.8 Tutorial)
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
## END FUNCTION ##

## Function to count the number of reads in which a primer is found
## (From DADA2 ITS Pipeline Workflow 1.8 Tutorial)
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE, max.mismatch = 2)
  return(sum(nhits > 0))
}
## END FUNCTION ##

## Function to count all primer orientations in a fastq file
## (From DADA2 ITS Pipeline Workflow 1.8 Tutorial)
## fn_r1: full filename of R1 read
## fn_r2: full filename of R2 read (complementary to R1; optional arg)
## primer_fwd: forward primer sequence
## primer_rev: reverse primer sequence
count_primer_orients <- function(fn_r1, fn_r2=NULL, primer_fwd, primer_rev) {
  fwd_orients <- allOrients(primer_fwd)
  rev_orients <- allOrients(primer_rev)
  if(is.null(fn_r2)) {
    rbind(FWDPrimer.R1.reads = sapply(fwd_orients, primerHits, fn = fn_r1), 
          REVPrimer.R1.reads = sapply(fwd_orients, primerHits, fn = fn_r1))
  } else {
    rbind(FWDPrimer.R1.reads = sapply(fwd_orients, primerHits, fn = fn_r1), 
          FWDPrimer.R2.reads = sapply(fwd_orients, primerHits, fn = fn_r2), 
          REVPrimer.R1.reads = sapply(rev_orients, primerHits, fn = fn_r1), 
          REVPrimer.R2.reads = sapply(rev_orients, primerHits, fn = fn_r2))
  }
}
## END FUNCTION ##

## Function to get sample names from list of file names, assuming sample name
## is the four underscore-separated pieces of the file name.
## Adapted from DADA2 ITS Tutorial (1.8)
get.sample.name <- function(fname) {
  paste(strsplit(basename(fname), "_")[[1]][1:4], collapse="_") # TODO: this could be made more robust
}
## END FUNCTION ##

# ## Function to get all unique filenames of samples with at least some
# ## reads remaining after the DADA2 pipeline.
# ## dir: Filenames from this directory should be listed
# ## trim_to_internal_lab_id: Whether to trim filenames to match format of
# ##      "internalLabID" in the mmg_soilRawDataFiles table
# get_fastq_names <- function(dir, trim_to_internal_lab_id=FALSE) {
#   f <- list.files(path = dir, pattern = ".fastq", full.names = TRUE)
#   if(trim_to_internal_lab_id) {
#     f <- basename(f)
#     f <- sub("^run[A-Za-z0-9]*_", "", f)
#     f <- sub("_((ITS)|(16S))_R[12].fastq(.tar)?(.gz)?(.tar)?", "", f)
#     return(f)
#   } else {
#     return(f)
#   }
# }
# ## END FUNCTION