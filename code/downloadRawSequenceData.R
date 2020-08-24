download_raw_sequence_files <- function(pathToMetadata, outDir = paste0(pathToMetadata, "/ECS_seqFiles"), 
    checkDownloadSize = TRUE, decompress = TRUE, keepZippedFiles = FALSE, removeTarFiles = TRUE) {
	# author: Lee F Stanish
	# Date: Aug 24, 2020
	# Description: download raw, per sample marker gene sequencing files from the mmg marker gene sequencing metadata
	# Requirements: mmg variables file, must be located in pathToMetadata/'variables.csv'
	# Requirements: mmg raw data files, must be located in pathToMetadata/ 
	# Output: Raw sequence files for every record listed in the metadata file provided in the pathToMetadata. If 
	# removeTarFiles=TRUE (recommended), any batch-level .tar.gz files will be removed from the 	
	# mmg_rawDataFiles table prior to downloading raw sequence data. This will overwrite the original mmg_rawDataFiles 
	# table. If removeTarFiles=FALSE, both per sample and batch-level data will be downloaded.
	# checkDownloadSize: indicates whether the size of the download is provided.
	# decompress: whether the downloaded data should be decompressed.
	# keepZippedFiles: indicates whether the originally downloaded, compressed files should be retained.
	
	library(neonUtilities)
	
	# check raw data files table for .tar.gz files and remove prior to downloading sequence data
	if(removeTarFiles==TRUE) {
		rawFilesList <- list.files(pathToMetadata, pattern="RawDataFiles.csv", full.names=TRUE)
		if(length(rawFilesList) > 1) {
			stop("More than 1 raw sequence metadata file found, only 1 file allowed")
		} else if(rawFilesList == "") {
			stop("No raw sequence metadata file found, at least 1 file is required")
		} else if(length(rawFilesList) == 1) {
			rawFilesDat <- read.csv(rawFilesList, stringsAsFactors=FALSE)
			tarGzInd <- grep('\\.tar\\.gz', rawFilesDat$rawDataFileName)
			if(length(tarGzInd) > 0) {
				rawFilesDat <- rawFilesDat[-tarGzInd, ]
				print(paste(length(tarGzInd), "row(s) containing batch-level sequence data were removed from 
				mmg_rawDataFiles dataset."))
				write.csv(rawFilesDat, rawFilesList, row.names=FALSE, na="")
			}
		}
		
	}
	zipsByURI(filepath=pathToMetadata, savepath = outDir, pick.files = FALSE, check.size = checkDownloadSize, unzip = 			decompress, saveZippedFiles = keepZippedFiles)
}

