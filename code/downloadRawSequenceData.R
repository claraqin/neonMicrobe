## Function downloads the metadata for various sequencing data products and downloads the ##
## raw data files ##

# Metadata is DP1.10108: "Soil microbe marker gene sequences"

# sites = "HARV"
# startYrMo = "2016-04"
# endYrMo = "2016-07"
# outdir = "/afs/cats.ucsc.edu/users/b/claraqin/zhulab/NEON_DoB_analysis/data/Illumina"
# # outdir = "~/Downloads"
# checkFileSize = TRUE

downloadRawSequenceData <- function(sites="all", startYrMo="YYYY-MM", endYrMo="YYYY-MM", 
                                    outdir="~/Downloads", checkFileSize=TRUE,
                                    return_data=TRUE) {
  # outdir - path to directory to ouput the data. Defaults to the R default directory if none provided
  # change checkFileSize to FALSE to override file size checks
  
  library(neonUtilities)
  library(utils)
  
  dat <- loadByProduct(dpID="DP1.10108.001", site=sites, startdate=startYrMo, enddate=endYrMo,
                       package="expanded", check.size=checkFileSize)
  
  u.urls <- unique(dat$mmg_soilRawDataFiles$rawDataFilePath)
  fileNms <- gsub('^.*\\/', "", u.urls)
  print(paste("There are", length(u.urls), " unique files to download.") )
  
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
  ## END FUNCTION ##
}
