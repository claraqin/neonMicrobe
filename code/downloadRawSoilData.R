## Function downloads NEON soil data to be associated with soil microbial
## community composition data

# Data product are
# - DP1.10078: "Soil chemical properties (Distributed periodic)"
# - DP1.10086: "Soil physical properties (Distributed periodic)"

# Used for testing:
# sites = "HARV"
# startYrMo = "2014-06"
# endYrMo = "2014-07"
# outdir = "/afs/cats.ucsc.edu/users/b/claraqin/zhulab/NEON_DoB_analysis/data/metadata"
# # outdir = "~/Downloads"
# checkFileSize = TRUE

downloadRawSoilData <- function(sites="all", startYrMo="YYYY-MM", endYrMo="YYYY-MM", 
                                    outdir="~/Downloads", checkFileSize=TRUE,
                                    return_data=TRUE) {
  # outdir - path to directory to ouput the data. Defaults to the R default directory if none provided
  # change checkFileSize to FALSE to override file size checks
  
  library(neonUtilities)
  library(dplyr)
  
  neonUtilities::zipsByProduct("DP1.10078.001", site=sites, startdate=startYrMo, 
                               enddate=endYrMo, check.size=TRUE, savepath = outdir)
  
  neonUtilities::zipsByProduct("DP1.10086.001", site=sites, startdate=startYrMo, 
                               enddate=endYrMo, check.size=TRUE, savepath = outdir)
  
  stackDir_chem <- paste(outdir, "filesToStack10078", sep="/")
  stackDir_phys <- paste(outdir, "filesToStack10086", sep="/")
  stackByTable(stackDir_chem, folder = TRUE)
  stackByTable(stackDir_phys, folder = TRUE)
  
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
    
    select(dat_soilCoreCollection, uid:biomassID) %>%
      left_join(select(dat_soilMoisture, uid:dryMassFraction)) %>%
      left_join(select(dat_soilpH, uid:caclpHRatio)) %>%
      left_join(select(dat_soilChemistry, uid:CNratio)) -> 
      dat_soil
    
    return(dat_soil)
  }
  ## END FUNCTION ##
}
