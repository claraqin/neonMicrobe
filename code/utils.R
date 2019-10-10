# Utilities for use in other R scripts in NEON/DoB microbial community analysis

library(neonUtilities)
library(utils)
library(dplyr)

# Preset parameters for the "All" variants of download functions:
preset_sites = "all"
preset_startYrMo = "2011-06"
preset_endYrMo = format(Sys.Date(), "%Y-%m")
preset_outdir_sequence = "/afs/cats.ucsc.edu/users/b/claraqin/zhulab/NEON_DoB_analysis/data/Illumina/NEON"
preset_outdir_seqmeta = "/afs/cats.ucsc.edu/users/b/claraqin/zhulab/NEON_DoB_analysis/data/sequence_metadata"
preset_outdir_soil = "/afs/cats.ucsc.edu/users/b/claraqin/zhulab/NEON_DoB_analysis/data/soil"
preset_checkFileSize = FALSE
preset_return_data = TRUE

## Function downloads the metadata for NEON marker gene sequencing data products 
#
# Metadata is just one table withing DP1.10108 ("Soil microbe marker gene sequences")
# called "mmg_soilRawDataFiles.csv"
downloadSequenceMetadata <- function(sites = preset_sites, startYrMo = preset_startYrMo, endYrMo = preset_endYrMo, 
                                     outdir = preset_outdir_seqmeta, checkFileSize = preset_checkFileSize, return_data = preset_return_data) {
  # outdir - path to directory to download the data
  # change checkFileSize to FALSE to override file size checks
  
  PRNUM <- "10108"
  dpID <- paste("DP1", PRNUM, "001", sep=".")
  stackDir <- paste(outdir, paste0("filesToStack",PRNUM), sep="/")
  
  # Download only if not already downloaded to outdir
  if(!dir.exists(stackDir)) {
    neonUtilities::zipsByProduct(dpID, site=sites, startdate=startYrMo, 
                                 enddate=endYrMo, package="expanded", check.size=checkFileSize, 
                                 savepath=outdir)
    
    stackByTable(stackDir, folder = TRUE)
  } else {
    warning("Data product ", PRNUM, " has already been downloaded to ",
            outdir, ".\nEnsure that sites / date range of downloaded data was correct,\nor else remove them and re-run your command.")
  }
  
  if(return_data) { # If not, then simply downloads the data to the outdir
    path <- paste(stackDir, "stackedFiles", 
                  "mmg_soilRawDataFiles.csv", sep="/")
    dat <- read.delim(path, sep=",", stringsAsFactors=FALSE)
    
    return(dat)
  }
}
## END FUNCTION ##

## Function downloads ALL metadata for NEON marker gene sequencing data products,
## using preset parameters
# 
# Wrapper for downloadSequenceMetadata,
# parameterized for socs-stats.ucsc.edu server,
# and will pull data up to current month.
downloadAllSequenceMetadata <- function(sites = preset_sites, startYrMo = preset_startYrMo, endYrMo = preset_endYrMo, 
                                        outdir = preset_outdir_seqmeta, checkFileSize = preset_checkFileSize, return_data = preset_return_data) {
  return(downloadSequenceMetadata(sites, startYrMo, endYrMo, 
                                  outdir, checkFileSize, return_data))
}

## Function downloads the metadata for NEON marker gene sequencing data products 
## AND downloads the NEON raw sequence data files
#
# Wrapper for downloadSequenceMetadata
downloadRawSequenceData <- function(sites = preset_sites, startYrMo = preset_startYrMo, endYrMo = preset_endYrMo, 
                                    outdir = preset_outdir_sequence, checkFileSize = preset_checkFileSize, return_data = preset_return_data) {
  # outdir - path to directory to download the data
  # change checkFileSize to FALSE to override file size checks
  
  metadata <- downloadSequenceMetadata(sites, startYrMo, endYrMo, outdir, checkFileSize=FALSE, return_data=TRUE)
  
  # u.urls <- unique(dat$mmg_soilRawDataFiles$rawDataFilePath)
  u.urls <- unique(metadata$rawDataFilePath)
  fileNms <- gsub('^.*\\/', "", u.urls)
  print(paste("There are", length(u.urls), "unique files to download.") )
  
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
    return(dat)
  }
}
## END FUNCTION ##


## Function downloads ALL metadata for NEON marker gene sequencing data products
## AND downloads ALL NEON raw sequence data, using preset parameters
#
# Wrapper for downloadRawSequenceData,
# parameterized for socs-stats.ucsc.edu server,
# and will pull data up to current month.
downloadAllRawSequenceData <- function(sites = preset_sites, startYrMo = preset_startYrMo, endYrMo = preset_endYrMo, 
                                       outdir = preset_outdir_sequence, checkFileSize = preset_checkFileSize, return_data = preset_return_data) {
  return(downloadRawSequenceData(sites, startYrMo, endYrMo, 
                                 outdir, checkFileSize, return_data))
}
## END FUNCTION


## Function downloads NEON soil data to be associated with soil microbial
## community composition data
#
# Data product are
# - DP1.10078: "Soil chemical properties (Distributed periodic)"
# - DP1.10086: "Soil physical properties (Distributed periodic)"
downloadRawSoilData <- function(sites = preset_sites, startYrMo = preset_startYrMo, endYrMo = preset_endYrMo, 
                                outdir = preset_outdir_soil, checkFileSize = preset_checkFileSize, return_data = preset_return_data) {
  # outdir - path to directory to download the data. Defaults to the R default directory if none provided
  # change checkFileSize to FALSE to override file size checks
  
  PRNUM_chem <- "10078"
  PRNUM_phys <- "10086"
  dpID_chem <- paste("DP1", PRNUM_chem, "001", sep=".")
  dpID_phys <- paste("DP1", PRNUM_phys, "001", sep=".")
  stackDir_chem <- paste(outdir, paste0("filesToStack",PRNUM_chem), sep="/")
  stackDir_phys <- paste(outdir, paste0("filesToStack",PRNUM_phys), sep="/")
  
  # Download only if not already downloaded to outdir
  if(!dir.exists(stackDir_chem)) {
    neonUtilities::zipsByProduct(dpID_chem, site=sites, startdate=startYrMo, 
                                 enddate=endYrMo, package="expanded", check.size=checkFileSize, savepath = outdir)
    
    stackByTable(stackDir_chem, folder = TRUE)
  } else {
    warning("Data product ", PRNUM_chem, " has already been downloaded to ",
            outdir, ".\nEnsure that sites / date range of downloaded data was correct,\nor else remove them and re-run your command.")
  }
  
  # Download only if not already downloaded to outdir
  if(!dir.exists(stackDir_phys)) {
    neonUtilities::zipsByProduct(dpID_phys, site=sites, startdate=startYrMo, 
                                 enddate=endYrMo, package="expanded", check.size=checkFileSize, savepath = outdir)
    
    stackByTable(stackDir_phys, folder = TRUE)
  } else {
    warning("Data product ", PRNUM_phys, " has already been downloaded to ",
            outdir, ".\nEnsure that sites / date range of downloaded data was correct,\nor else remove them and re-run your command.")
  }
  
  if(return_data) { # If not, then simply downloads the data to the outdir
    
    path_soilCoreCollection <- paste(stackDir_phys, "stackedFiles", 
                                     "sls_soilCoreCollection.csv", sep="/")
    path_soilMoisture <- paste(stackDir_phys, "stackedFiles", 
                               "sls_soilMoisture.csv", sep="/")
    path_soilpH <- paste(stackDir_phys, "stackedFiles", 
                         "sls_soilpH.csv", sep="/")
    # path_bgcSubsampling <- paste(stackDir_phys, "stackedFiles", 
    #                              "sls_bgcSubsampling.csv", sep="/")
    path_soilChemistry <- paste(stackDir_chem, "stackedFiles", 
                                "sls_soilChemistry.csv", sep="/")
    
    dat_soilCoreCollection <- read.delim(path_soilCoreCollection, sep=",", stringsAsFactors=FALSE)
    dat_soilMoisture <- read.delim(path_soilMoisture, sep=",", stringsAsFactors=FALSE)
    dat_soilpH <- read.delim(path_soilpH, sep=",", stringsAsFactors=FALSE)
    dat_soilChemistry <- read.delim(path_soilChemistry, sep=",", stringsAsFactors=FALSE)
    
    joining_cols <- c("domainID", "siteID", "plotID", "sampleID")
    
    select(dat_soilCoreCollection, uid:biomassID, -uid) %>% 
      full_join(select(dat_soilMoisture, uid:dryMassFraction, -uid), by=joining_cols) %>%
      full_join(select(dat_soilpH, uid:caclpHRatio, -uid), by=joining_cols) %>%
      full_join(select(dat_soilChemistry, uid:CNratio, -uid), by=joining_cols) -> 
      dat_soil
    
    return(dat_soil)
  }
}
## END FUNCTION ##


## Function downloads ALL NEON soil data, using preset parameters
#
# Wrapper for downloadRawSoilData,
# parameterized for socs-stats.ucsc.edu server,
# and will pull data up to current month.
downloadAllRawSoilData <- function(sites = preset_sites, startYrMo = preset_startYrMo, endYrMo = preset_endYrMo, 
                                   outdir = preset_outdir_soil, checkFileSize = preset_checkFileSize, return_data = preset_return_data) {
return(downloadRawSoilData(sites, startYrMo, endYrMo, 
                             outdir, checkFileSize, return_data))
}
## END FUNCTION ##
