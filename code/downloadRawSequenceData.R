download_raw_sequence_files <- function(meta="", pathToMetadata=file.path(PRESET_OUTDIR, PRESET_OUTDIR_SEQMETA), 
                                        outDir = file.path(PRESET_OUTDIR, PRESET_OUTDIR_SEQUENCE), checkDownloadSize = TRUE, 
                                        decompress = TRUE, keepZippedFiles = FALSE, removeTarFiles = TRUE) {
	# author: Lee F Stanish
	# Date: Dec 4, 2020
	# Description: download raw, per sample marker gene sequencing files from the mmg marker gene sequencing metadata
	# Requirements: mmg variables file, must be located in pathToMetadata/'variables.csv'
	# Requirements: mmg raw data files, must be located in pathToMetadata/ 
	# Output: Raw sequence files for every record listed in the metadata file provided in the pathToMetadata. If 
	# removeTarFiles=TRUE (recommended), any batch-level .tar.gz files will be removed from the 	
	# mmg_rawDataFiles table prior to downloading raw sequence data. This will overwrite the original mmg_rawDataFiles 
	# table. If removeTarFiles=FALSE, both per sample and batch-level data will be downloaded.
	# checkDownloadSize: indicates whether the size of the download is provided.
	# decompress: whether the downloaded data should be decompressed (default=TRUE).
	# keepZippedFiles: indicates whether the originally downloaded, compressed files should be retained (default=FALSE).
	
	library(neonUtilities)
	
  metaFilesList <- list.files(pathToMetadata, pattern="Metadata|metadata", full.names=TRUE)
  
	# define dataset
  if(meta != "") {
    dat <- meta
  } else {
    if(length(metaFilesList) > 1) {
      stop("More than 1 raw sequence metadata file found, only 1 file allowed")
    } else if(metaFilesList == "") {
      stop("No raw sequence metadata file found, at least 1 file is required")
    } else if(length(metaFilesList) == 1) {
      dat <- read.csv(metaFilesList, stringsAsFactors=FALSE)
    }
  }
  
  # Remove tarfiles if specified, export modified metadata table
  tarGzInd <- grep('\\.tar\\.gz', dat$rawDataFileName)
  if(length(tarGzInd) > 0) {
    dat <- dat[-tarGzInd, ]
    write.csv(dat, paste0(pathToMetadata, "/metadata_batchFilesRemoved.csv"), row.names=FALSE, na="")
    print(paste(length(tarGzInd), "row(s) containing batch-level sequence data removed from dataset. Updated metadata file downloaded to:",
                paste0(pathToMetadata, "/metadata_batchFilesRemoved.csv")) )
  }
      
	zipsByURI(filepath=pathToMetadata, savepath = outDir, pick.files = FALSE, check.size = checkDownloadSize, unzip = decompress, saveZippedFiles = keepZippedFiles)
}

