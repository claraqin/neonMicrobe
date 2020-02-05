# Utilities for use in other R scripts in NEON/DoB microbial community analysis

library(neonUtilities)
library(utils)
library(dplyr)

# Load parameters from params.R
source("./code/params.R")

# If preset output directories do not exist, create them
if(!dir.exists(preset_outdir_sequence)) dir.create(preset_outdir_sequence, recursive=TRUE)
if(!dir.exists(preset_outdir_seqmeta)) dir.create(preset_outdir_seqmeta, recursive=TRUE)
if(!dir.exists(preset_outdir_soil)) dir.create(preset_outdir_soil, recursive=TRUE)

## Function to write site and date range parameters. To be called within other functions.
#
write_site_and_date_range <- function(writeDir, sites, startYrMo, endYrMo) {
  # paste("writing in ", writeDir)
  write.table(
    data.frame(x1=c("sites", "startYrMo", "endYrMo"),
               x2=c(sites, startYrMo, endYrMo)),
    file = paste(writeDir, site_and_date_range_filename, sep="/"),
    sep=":", col.names=FALSE, quote=FALSE, row.names=FALSE
  )
}
## END FUNCTION


## Function issues warning that metadata has already been downloaded, and indicates the sites
## and date ranges for which it was downloaded. To be called within other functions.
#
warn_already_downloaded <- function(PRNUM, outdir) {
  # print(paste("checking in ", outdir))
  stackDir <- paste(outdir, paste0("filesToStack",PRNUM), sep="/")
  # print(paste("specifically in ", stackDir))
  site_and_date_range <- read.table(paste(stackDir, site_and_date_range_filename, sep="/"), 
                                    header=FALSE, sep=":")
  warning("Data product ", PRNUM, 
          " has already been downloaded to ",outdir, 
          ".\nEnsure that sites / date range of downloaded data was correct,\n",
          "or else remove the data and re-run your command.\n",
          "sites:", site_and_date_range[1,2],
          "  startYrMo:", site_and_date_range[2,2],
          "  endYrMo:", site_and_date_range[3,2])
}
## END FUNCTION


## Function downloads the metadata for NEON marker gene sequencing data products 
#
# Metadata is just one table withing DP1.10108 ("Soil microbe marker gene sequences")
# called "mmg_soilRawDataFiles.csv"
downloadSequenceMetadata <- function(sites = preset_site, startYrMo = preset_startYrMo, endYrMo = preset_endYrMo, 
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
    
    print(paste("Attempt to stackByTable in", stackDir))
    stackByTable(stackDir, folder = TRUE)
    
    # Write site_and_date_range.txt file to record parameters
    write_site_and_date_range(stackDir, sites, startYrMo, endYrMo)
    
  } else {
    warn_already_downloaded(PRNUM, outdir)
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
downloadAllSequenceMetadata <- function() {
  return(downloadSequenceMetadata(preset_site, preset_startYrMo, preset_endYrMo, 
                                  preset_outdir_seqmeta, preset_checkFileSize, preset_return_data))
}

## Function downloads the metadata for NEON marker gene sequencing data products 
## AND downloads the NEON raw sequence data files
#
# Wrapper for downloadSequenceMetadata
downloadRawSequenceData <- function(sites = preset_site, startYrMo = preset_startYrMo, endYrMo = preset_endYrMo, 
                                    outdir = preset_outdir_sequence, checkFileSize = preset_checkFileSize, return_data = preset_return_data) {
  # outdir - path to directory to download the data
  # change checkFileSize to FALSE to override file size checks
  
  metadata <- downloadSequenceMetadata(sites, startYrMo, endYrMo, checkFileSize=FALSE, return_data=TRUE)
  
  u.urls <- unique(metadata$rawDataFilePath)
  fileNms <- gsub('^.*\\/', "", u.urls)
  print(paste("There are", length(u.urls), "unique raw sequence files to download.") )
  
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


## Function downloads ALL metadata for NEON marker gene sequencing data products
## AND downloads ALL NEON raw sequence data, using preset parameters
#
# Wrapper for downloadRawSequenceData,
# parameterized for socs-stats.ucsc.edu server,
# and will pull data up to current month.
downloadAllRawSequenceData <- function() {
  return(downloadRawSequenceData(preset_site, preset_startYrMo, preset_endYrMo, 
                                 preset_outdir_sequence, preset_checkFileSize, preset_return_data))
}
## END FUNCTION


## Function downloads NEON soil data to be associated with soil microbial
## community composition data
#
# Data product are
# - DP1.10078: "Soil chemical properties (Distributed periodic)"
# - DP1.10086: "Soil physical properties (Distributed periodic)"
downloadRawSoilData <- function(sites = preset_site, startYrMo = preset_startYrMo, endYrMo = preset_endYrMo, 
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
    
    # Write site_and_date_range.txt file to record parameters
    write_site_and_date_range(stackDir_chem, sites, startYrMo, endYrMo)
    
  } else {
    warn_already_downloaded(PRNUM_chem, outdir)
  }
  
  # Download only if not already downloaded to outdir
  if(!dir.exists(stackDir_phys)) {
    neonUtilities::zipsByProduct(dpID_phys, site=sites, startdate=startYrMo, 
                                 enddate=endYrMo, package="expanded", check.size=checkFileSize, savepath = outdir)
    
    stackByTable(stackDir_phys, folder = TRUE)
    
    # Write site_and_date_range.txt file to record parameters
    write_site_and_date_range(stackDir_phys, sites, startYrMo, endYrMo)
    
  } else {
    warn_already_downloaded(PRNUM_phys, outdir)
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
downloadAllRawSoilData <- function() {
return(downloadRawSoilData(preset_site, preset_startYrMo, preset_endYrMo, 
                           preset_outdir_soil, preset_checkFileSize, preset_return_data))
}
## END FUNCTION ##
