## Function downloads the metadata for various sequencing data products and downloads the ##
## raw data files ##

downloadRawSequenceData <- function(sites="all", startYrMo="YYYY-MM", endYrMo="YYYY-MM", 
                                    outdir="~/Downloads", checkFileSize=TRUE) {
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
                                              paste(outdir, fileNms, sep="/"), 
                                             paste(getwd(), fileNms, sep="/" )) ) 
    if(dir.exists(outdir)) {
      print(paste("Finished downloading", paste(outdir, fileNms, sep="/")))
    } else {
      print(paste("Finished downloading", paste(getwd(), fileNms, sep="/" )))
    }
  }
  
  return(dat)
  ## END FUNCTION ##
}
