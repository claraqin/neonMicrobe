# Utility functions for the neonMicrobe package

# This version of "utils.R" is no longer the only script containing
# functions from this package. The functions have now been spread
# across multiple R scripts in the "R" subdirectory. The functions
# remaining in this script are really "utilities" in the sense that
# they have more general purposes than the others in this package.

#' Remove Unmatched Fastq Files (DEPRECATED)
#'
#' NOTE: This function is now deprecated in favor of \code{\link{getPairedFastqFiles}} and \code{\link{removeUnpairedFastqFiles}}.
#'
#' Removes any forward-read files that do not have reverse-read counterparts, and vise versa.
#' This function is necessary because dada2::filterAndTrim() will throw an error
#' if the forward-read files and the reverse-read files are mismatched.
#'
#' @param fnFs Full name(s) of fastq file(s) containing forward-read sequence data.
#' @param fnRs Full name(s) of fastq file(s) containing reverse-read sequence data.
#' @param post_samplename_pattern (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#' @param verbose Whether to print messages regarding which files are matched and which are unmatched.
#'
#' @return List of length 2. The first element is a vector of forward-read files that have reverse-read counterparts; the second element is a vector of reverse-read files that have forward-read counterparts.
#' @export
#'
#' @examples
#' matched_fn <- remove_unmatched_files(c("sample1_R1.fastq", "sample2_R1.fastq"), c("sample1_R2.fastq", "sample2_R2.fastq", "sample3_R2.fastq"))
#' fnFs <- matched_fn[[1]]
#' fnRs <- matched_fn[[2]]
remove_unmatched_files <- function(fnFs, fnRs, post_samplename_pattern = "_R(1|2).*\\.fastq", verbose=FALSE){
  message("remove_unmatched_files() is deprecated. Please use getPairedFastqFiles() or removeUnpairedFastqFiles() instead.")
  basefilenames_Fs <- sapply(strsplit(fnFs, post_samplename_pattern), `[`, 1)
  basefilenames_Rs <- sapply(strsplit(fnRs, post_samplename_pattern), `[`, 1)
  rm_from_fnFs <- which(!(basefilenames_Fs %in% basefilenames_Rs))
  rm_from_fnRs <- which(!(basefilenames_Rs %in% basefilenames_Fs))

  if (length(c(rm_from_fnFs, rm_from_fnRs)) == 0){
    if(verbose == TRUE) message("All R1/R2 files had counterparts.\n")
  } else {
    for(i in rm_from_fnFs) {
      if(verbose) message(paste(basefilenames_Fs[i], "does not have an R2 counterpart. Omitting from this analysis.\n"))
      fnFs <- fnFs[-rm_from_fnFs]
    }
    for(i in rm_from_fnRs) {
      if(verbose) message(paste(basefilenames_Rs[i], "does not have an R1 counterpart. Omitting from this analysis.\n"))
      fnRs <- fnRs[-rm_from_fnRs]
    }
  }
  return(list(R1=fnFs, R2=fnRs))
}


#' Get Paired Fastq Files
#'
#' Given a vector of fastq file names and corresponding metadata, outputs
#' a list containing matching R1 and R2 files in the same order.
#' Useful for providing files to functions that utilize
#' \code{\link[dada2]{filterAndTrim}}, such as \code{\link{trimPrimer16S}}
#' and \code{\link{runDadaITS}}.
#'
#' See related function: \code{\link{removeUnpairedFastqFiles}}.
#'
#' @param fn Full name(s) of fastq file(s). Tolerant to filenames with appended run ID prefix (e.g. "runB69RN_...") and agnostic to .gz compression.
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param value (Optional) Default TRUE; returns file names. If FALSE, returns indices.
#' @param verbose (Optional) Default TRUE. Whether to print warning messages when metadata does not contain records corresponding to all provided fastq files.
#'
#' @return List of length 2; first element contains matched R1 (forward) file names, and second element contains corresponding R2 (reverse) file names.
#' @export
#'
#' @examples
#' \dontrun{
#' #' matched_fn <- getPairedFastqFiles(c("sample1_R1.fastq", "sample3_R1.fastq", "sample1_R2.fastq", "sample2_R2.fastq", "sample3_R2.fastq"), meta)
#' fnFs <- matched_fn[[1]] # "sample1_R1.fastq" "sample3_R1.fastq"
#' fnRs <- matched_fn[[2]] # "sample1_R2.fastq" "sample3_R2.fastq"
#' }
#' @importFrom magrittr "%>%"
getPairedFastqFiles <- function(fn, meta, value=TRUE, verbose=TRUE) {
  meta_ext <- matchFastqToMetadata(fn, meta, verbose=verbose)
  meta_ext$orientation <- if_else(grepl("R1", meta_ext$rawDataFileDescription), "R1",
                                  if_else(grepl("R2", meta_ext$rawDataFileDescription), "R2",
                                          NA_character_))

  meta_ext %>%
    dplyr::filter(!is.na(uid.rawFiles)) %>%
    dplyr::group_by(sequencerRunID, dnaSampleID) %>%
    dplyr::summarise(n_orientations = n_distinct(orientation, na.rm=TRUE), .groups="drop") ->
    meta_summ
  if(!any(meta_summ$n_orientations == 2)) {
    message("There were no matches. Check your inputs and try again.")
    return(NULL)
  } else {
    matched_ids <- meta_summ$dnaSampleID[meta_summ$n_orientations==2]
  }

  if(value==TRUE) {
    return(list(
      meta_ext$file[meta_ext$orientation=="R1" & meta_ext$dnaSampleID %in% matched_ids],
      meta_ext$file[meta_ext$orientation=="R2" & meta_ext$dnaSampleID %in% matched_ids]
    ))
  } else {
    return(list(
      which(meta_ext$orientation=="R1" & meta_ext$dnaSampleID %in% matched_ids),
      which(meta_ext$orientation=="R2" & meta_ext$dnaSampleID %in% matched_ids)
    ))
  }
}


#' Remove Unpaired Fastq Files
#'
#' Given a vector of fastq file names and corresponding metadata, outputs
#' a list containing matching R1 and R2 files in the same order.
#' Useful for providing files to functions that utilize
#' \code{\link[dada2]{filterAndTrim}}, such as \code{\link{trimPrimer16S}}
#' and \code{\link{runDadaITS}}.
#'
#' See related function: \code{\link{getPairedFastqFiles}}.
#'
#' @param fnFs,fnRs Full name(s) of fastq file(s) corresponding to R1 (forward) reads and R2 (reverse) reads, respectively. Tolerant to filenames with appended run ID prefix (e.g. "runB69RN_...") and agnostic to .gz compression.
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param value (Optional) Default TRUE; returns file names. If FALSE, returns indices.
#' @param verbose (Optional) Default TRUE. Whether to print warning messages when metadata does not contain records corresponding to all provided fastq files.
#'
#' @return List of length 2; first element contains matched R1 (forward) file names, and second element contains corresponding R2 (reverse) file names.
#' @export
#'
#' @examples
#' \dontrun{
#' #' matched_fn <- removeUnpairedFastqFiles(c("sample1_R1.fastq", "sample3_R1.fastq"), c("sample1_R2.fastq", "sample2_R2.fastq", "sample3_R2.fastq"), meta)
#' fnFs <- matched_fn[[1]] # "sample1_R1.fastq" "sample3_R1.fastq"
#' fnRs <- matched_fn[[2]] # "sample1_R2.fastq" "sample3_R2.fastq"
#' }
#'
#' @importFrom magrittr "%>%"
removeUnpairedFastqFiles <- function(fnFs, fnRs, meta, value=TRUE, verbose=TRUE) {
  meta_ext <- matchFastqToMetadata(c(fnFs, fnRs), meta, verbose=verbose)
  meta_ext$orientation <- c(rep("R1", length(fnFs)), rep("R2", length(fnRs)))

  meta_ext %>%
    dplyr::filter(!is.na(uid.rawFiles)) %>%
    dplyr::group_by(sequencerRunID, dnaSampleID) %>%
    dplyr::summarise(n_orientations = n_distinct(orientation, na.rm=TRUE), .groups="drop") -> meta_summ
  if(!any(meta_summ$n_orientations == 2)) {
    message("There were no matches. Check your inputs and try again.")
    return(NULL)
  } else {
    matched_ids <- meta_summ$dnaSampleID[meta_summ$n_orientations==2]
  }

  if(value==TRUE) {
    return(list(
      meta_ext$file[meta_ext$orientation=="R1" & meta_ext$dnaSampleID %in% matched_ids],
      meta_ext$file[meta_ext$orientation=="R2" & meta_ext$dnaSampleID %in% matched_ids]
    ))
  } else {
    return(list(
      which(meta_ext$orientation=="R1" & meta_ext$dnaSampleID %in% matched_ids),
      which(meta_ext$orientation=="R2" & meta_ext$dnaSampleID %in% matched_ids)
    ))
  }
}


#' Match Fastq Files to Metadata
#'
#' Helper function to match fastq file names to sequence metadata via
#' the sequence metadata table's "rawDataFileName" column. Tolerant to
#' filenames with appended run ID prefix (e.g. "runB69RN_...") and
#' agnostic to .gz compression.
#'
#' Used in \code{\link{getPairedFastqFiles}} and \code{\link{removeUnpairedFastqFiles}}.
#'
#' @param fn Full name(s) of fastq file(s).
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param verbose (Optional) Default TRUE. Whether to print warning messages when metadata does not contain records corresponding to all provided fastq files.
#'
#' @return Data frame, in which the first column contains the fastq file names in the order received, and the remaining columns contain corresponding metadata records. Has as many rows as the number of input filenames.
#' @export
#'
#' @examples
#' \dontrun{
#' meta_ext <- matchFastqToMetadata(c("sample1_R1.fastq", "sample1_R2.fastq"), meta)
#' }
matchFastqToMetadata <- function(fn, meta, verbose=TRUE) {
  # Get basenames of fn
  fn_base <- basename(fn)

  # Remove runID if appended to beginning of filename
  key <- sub("^run[A-Za-z0-9]*_", "", fn_base)

  # Append ".gz" to end of filename if missing
  key[!grepl(".gz$", key)] <- paste0(key[!grepl(".gz$", key)], ".gz")

  key_match <- match(key, as.character(meta$rawDataFileName))
  if(any(is.na(key_match))) {
    if(verbose) {
      message("Matching file names to metadata: ", sum(is.na(key_match)), " files did not have matching records in the provided metadata. ",
              "Double-check to ensure that the provided metadata is of the appropriate scope.")
    }
  }
  return(cbind(file = fn, meta[key_match,], stringsAsFactors=FALSE))
}


#' Get Truncation Length
#'
#' Decides on truncation length for trimmed 16S reads based on quality score means. Default cutoff is a score of 30. Warns you if the beginning (first 10) bases are low-quality, but returns the first low-quality base after the 10th base. If no bases are below score, returns the last base. The truncation lengths can be aggregated (e.g. minimum) and used as an argument in \code{\link[dada2]{filterAndTrim}}.
#'
#' @param fl Full names of fastq files.
#' @param qscore Default 30. Mean quality score threshold at which to truncate the remainder of the read.
#' @param n Default 5e+05. The number of reads to sample when processing fastq files.
#' @param verbose Default TRUE. Whether to return message regarding truncation length for each file. Includes warning messages.
#'
#' @return Integer vector of truncation lengths to use for the set of fastq files.
#' @export
#'
#' @examples
#' \dontrun{
#' trunc_len_fwd <- getTruncationLength(c("sample1_R1.fastq", "sample2_R1.fastq"))
#' trunc_len_rev <- getTruncationLength(c("sample1_R2.fastq", "sample2_R2.fastq"))
#' }
getTruncationLength <-function (fl, qscore = 30, n = 5e+05, verbose = TRUE){
  trunc_lengths <- data.frame(file = character(0), early_lowqual = numeric(0),
                              trunc_length = numeric(0))
  for (f in fl) {
    srqa <- ShortRead::qa(f, n = n)
    df <- srqa[["perCycle"]]$quality
    means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count,
                                                          df$Cycle)
    lowqual <- which(means < qscore)
    n_lowqual_in_first_20 <- sum(lowqual <= 20)
    first_lowqual_after_20 <- lowqual[lowqual > 20][1]

    trunc_lengths <- rbind(trunc_lengths, data.frame(file = f,
                                                     early_lowqual = n_lowqual_in_first_20,
                                                     trunc_length = ifelse(is.na(first_lowqual_after_20), length(means), first_lowqual_after_20)))

    if (verbose == TRUE) {
      if(lowqual_in_first_20 > 0){
        cat(paste0("Note: for sample ", basename(f),': ', n_lowqual_in_first_20, ' base(s) within first 20 bases are below your quality score. Consider trimming left.\n'))
      }
      if (is.na(lowqual_after_20)){
        cat(paste0(basename(f), ': After first 20 bases, no bases with a mean under your quality score found. Truncate at end of read, base: ', length(means),'\n'))

      } else if (!is.na(lowqual_after_20)){
        cat(paste0(basename(f),': After first 20 bases, first mean below your quality score is base: ',first_lowqual_after_20,'\n'))

      } else "Something's wrong. Inspect function/reads."
    } # end printing to console
  } # end loop
  return(trunc_lengths$trunc_length)
}



#' Get Number of Reads
#'
#' Gets the total number of reads in a dada-class object.
#' Primarily used for tracking reads throughout the pipeline.
#'
#' @param x A dada-class object.
#'
#' @return Integer; total number of reads in the dada-class object.
getN <- function(x) {
  sum(getUniques(x))
}


#' Get Sample Name from Fastq Filename (DEPRECATED)
#'
#' Deprecated in favor of code that reads dnaSampleID from sequence metadata.
#'
#' Extracts the sample names from fastq filenames, given regular expressions
#' denoting the different endings of R1 and R2 filenames.
#'
#' @param fn Fastq filename(s).
#' @param post_samplename_pattern1,post_samplename_pattern2 Default "_R1.*\\.fastq" and "_R2.*\\.fastq". Regular expressions denoting the ends of the R1 and R2 filenames, respectively, which will be trimmed away, leaving identical names between corresponding R1 and R2 files.
#'
#' @return Character vector; sample names of the input fastq filenames.
#' @export
getSampleName <- function(fn, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq") {
  message("getSampleName() is deprecated. Please retrieve dnaSampleID from sequence metadata instead.")
  sapply(strsplit(basename(fn), paste0("(",post_samplename_pattern1,")|(",post_samplename_pattern2, ")")), `[`, 1)
}


#' Parse Parameter Values from Rownames
#'
#' Parses parameter values from rownames. Particularly useful for sensitivity analyses, when rows have been constructed to contain information about the parameter values, e.g. "maxEE.R_8_truncLen.R_250".
#'
#' @param df Data frame with rownames containing parameter values of the form [param_name1]_[param_value1]_[param_name2]_[param_value2].
#' @param keep_rownames (Default FALSE) Whether to keep the rownames of the original dataframe.
#' @param keep_original_cols (Default TRUE) Whether to keep the columns of the original dataframe.
#' @param param1,param2 Parameter names embedded within the rownames of df.
#'
#' @return A data frame containing columns corresponding to the values of the parameters of interest, extracted from the rownames of the original data frame.
#' @export
#'
#' @examples \dontrun{ trackReadsWithParamVals <- parseParamsFromRownames(trackReads, PARAM1, PARAM2) }
parseParamsFromRownames <- function(df, param1, param2, keep_rownames=FALSE, keep_original_cols=TRUE) {
  if(class(df)=="matrix") df <- as.data.frame(df)
  param1.0 <- as.numeric(unlist(regmatches(rownames(df), gregexpr(paste0(param1, "_\\K[0-9]*"), rownames(df), perl=TRUE))))
  param2.0 <- as.numeric(unlist(regmatches(rownames(df), gregexpr(paste0(param2, "_\\K[0-9]*"), rownames(df), perl=TRUE))))
  runID0 <- unlist(regmatches(rownames(df), gregexpr("\\|\\Krun[A-Za-z0-9]*", rownames(df), perl=TRUE)))
  sampleID0 <- unlist(regmatches(rownames(df), gregexpr("\\|\\K.*", rownames(df), perl=TRUE)))
  if(length(param1.0)==nrow(df)) df[[param1]] <- param1.0
  if(length(param2.0)==nrow(df)) df[[param2]] <- param2.0
  if(length(runID0)==nrow(df)) df[["runID"]] <- runID0
  if(length(sampleID0)==nrow(df)) df[["sampleID"]] <- sampleID0
  if(!keep_rownames) {
    rownames(df) <- NULL
  }
  if(keep_original_cols) {
    return(df)
  } else {
    param1 <- enquo(param1)
    param2 <- enquo(param2)
    return(dplyr::select(df, any_of(c(!!param1, !!param2, "runID", "sampleID"))))
  }
}


#' Read Sequence Metadata into Function
#'
#' Helper function to accept two formats of metadata into various functions in this package: data.frame and filepath to a csv
#'
#' @param metadata Either a data.frame returned from \code{\link{downloadSequenceMetadata}} or the filepath to a local csv copy of the output from \code{\link{downloadSequenceMetadata}}.
#'
#' @return data.frame of the metadata.
readSequenceMetadata <- function(metadata) {
  metadata_load_err <- FALSE
  if(class(metadata) == "data.frame") return(metadata)
  if(class(metadata) == "character") {
    if(file.exists(metadata)) {
      return(read.csv(metadata))
    }
  }
  stop("'metadata' must be the data.frame output from downloadSequenceMetadata() or ",
       "the filepath to a local csv copy of the output from downloadSequenceMetadata()")
}

#' Combine 16S Read-Tracking Tables
#'
#' Convenience function for combining the read-tracking tables
#' across the trimming, filtering, and DADA steps of the 16S
#' processing pipeline.
#'
#' @param trim_table Read-tracking table output by \code{\link{trimPrimers16S}}.
#' @param filter_table Read-tracking table output by \code{\link{qualityFilter16S}}.
#' @param dada_table Read-tracking table output by \code{\link{runDada16S}}.
#'
#' @return An integer matrix with columns representing each step in the DADA2
#' 16S processing pipeline, and rows representing each input sample,
#' with values representing the read counts remaining at that step
#' for each sample.
#' @export
combineReadTrackingTables16S <- function(trim_table, filter_table, dada_table) {
  track <- Reduce(
    function(x, y, ...) transform(merge(x, y, by = 0, all = TRUE, ...), row.names=Row.names, Row.names = NULL),
    list(trim_table,
         filter_table[,2,drop=FALSE],
         dada_table)
  )
  names(track)[3] <- "filtered"
  track[is.na(track)] <- 0

  # # optionally save the combined table
  # if(!identical(out_file, FALSE)) {
  #   if(is.null(out_file)) {
  #     out_file <- file.path(NEONMICROBE_DIR_TRACKREADS(), "16S",
  #                           paste0("track_16s_", sub(" ", "_", gsub(":", "", Sys.time())), ".csv"))
  #   }
  #   write.csv(track, file = out_file)
  #   message("Saved read-tracking table to ", out_file)
  # }
  return(track)
}

#' Combine ITS Read-Tracking Tables
#'
#' Convenience function for combining the read-tracking tables
#' across the trimming, filtering, and DADA steps of the ITS
#' processing pipeline.
#'
#' @param trim_table Read-tracking table output by \code{\link{trimPrimersITS}}.
#' @param filter_table Read-tracking table output by \code{\link{qualityFilterITS}}.
#' @param dada_table Read-tracking table output by \code{\link{runDadaITS}}.
#'
#' @return An integer matrix with columns representing each step in the DADA2
#' ITS processing pipeline, and rows representing each input sample,
#' with values representing the read counts remaining at that step
#' for each sample.
#' @export
combineReadTrackingTablesITS <- function(prefilter_table, filter_table, dada_table) {
  track <- Reduce(
    function(x, y, ...) transform(merge(x, y, by=0, all = TRUE, ...), row.names=Row.names, Row.names=NULL),
    list(prefilter_table,
         filter_table,
         dada_table)
  )
  names(track)[1:4] <- c("reads.in", "prefiltered", "trimmed", "filtered")
  track[is.na(track)] <- 0

  return(track)
}
