#' QC Sequence Data
#'
#' Performs basic QAQC checks on sequence metadata prior to downloading sequence data and performing bioinformatics processing.
#' Running this function will remove metadata records for samples that do not meet user specifications. This will reduce the number of sequence files that are downloaded to only those that will be used for analysis, thereby saving file space and reducing download times.
#'
#' @param metadata The output of downloadSequenceMetadata(). Must be provided as either the data.frame returned by downloadSequenceMetadata() or as a filepath to the csv file produced by downloadSequenceMetadata() when outdir is provided.
#' @param outDir Default file.path(PRESET_OUTDIR_SEQMETA). Directory where raw sequence files can be found before reorganizing.
#' @param keepR2 "Y" (default) or "N". Should the reverse reads for a sample be retained? If "Y", then only samples that have both the forward (R1) and reverse (R2) reads will be retained.
#' @param rmDupes TRUE (default) or FALSE. Should duplicate samples be removed? If TRUE, then only the first records encountered for a particular sample identifier will be retained.
#'
#' @return Character vector of the files (including files within tarballs) that were successfully reorganized. If no files were successfully reorganized, returns no value.
qcMetadata <- function(metadata, outDir=getwd(), keepR2="Y", rmDupes=TRUE) {
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
  
  # validate keepR2
  if(!(keepR2 %in% c("Y", "N")) ) {
    stop("value for argument keepR2 invalid. Must be 'Y' or 'N.")
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
  
  # check for and remove duplicate sequence file names
  print("QC check for duplicate sequence file names...")
  # Add pause #
  Sys.sleep(1)
  dupeSeqIDs <- as.character(metadata$rawDataFileName[duplicated(metadata$rawDataFileName)] )
  if(length(dupeSeqIDs)==0) {
    print("QC check Pass. No duplicate sequence file names.")
  } else {
    numDupes <- length(dupeSeqIDs)
    print(paste0("QC check Fail. ", numDupes, " duplicate sequence file names found. Removing duplicated file(s). NOTE: files may not be identical. Review output before proceeding.") )
    print(paste0("Removing duplicated row: ", which(duplicated(metadata$rawDataFileName)) ) )
    metadata <- metadata[!duplicated(metadata$rawDataFileName), ]
  }
  
   
  # check for and flag duplicate dnaSampleIDs
  print("QC checking duplicate dnaSampleIDs...")
  # Add pause #
  Sys.sleep(1)
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
  Sys.sleep(1)
  # Handle sequence data with missing R2 data
  if(length(missingR2)>0) {
    print(paste0("Reverse read missing from ", length(missingR2), " records") )
    # If specified, remove R1 file
    if(keepR2=="Y") {
      print("Removing R1 files lacking a matching R2 file (default action when keepR2='Y')" )
      metaNotFlagged <- metaNotFlagged[-intersect(which(metaNotFlagged$runDir=="R1"), which(metaNotFlagged$dnaSampleID %in% missingR2) ), ]
    } else {
      metaNotFlagged$runDirFlag[intersect(which(metaNotFlagged$runDir=="R1"), which(metaNotFlagged$dnaSampleID %in% missingR2) )] <- "1"
      print("Flagging R1 files lacking a matching R2 file (default action when keepR2='N')" )
    }
  }
  # Handle sequence data with missing R1 data
  if(length(missingR1)>0) {
    print(paste0("Forward read missing from ", length(missingR1), " record(s). Removing R2 file for affected sequence data set(s).") )
    # remove R2 file
    metaNotFlagged <- metaNotFlagged[-intersect(which(metaNotFlagged$runDir=="R2"), which(metaNotFlagged$dnaSampleID %in% missingR1) ), ]
  }
  
  # Recombine original flagged records and remaining records post-initial flagging.
  out <- plyr::join(metaNotFlagged, metaFlagged)
  
  write.csv(out, paste0(qcDir, "mcc_metadata_", targetGene, "_QCed_", Sys.Date(), '.csv'), row.names = FALSE)
  return(print(paste("Output QCed file contains", nrow(out), "rows. File saved to the following directory:", qcDir)))
}
  
  
  
  
  
  
  
