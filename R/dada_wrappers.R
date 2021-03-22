#' Trim Primers from 16S Sequences
#'
#' Trims primers from 16S sequences using \code{\link[dada2]{filterAndTrim}}. This function
#' assumes that each read begins with its full primer sequence and operates
#' by truncating the beginning of each read by the length of its primer
#' sequence.
#'
#' @param fn Names of input fastq files, excluding directory path which is specified by dir_in. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param dir_in Directory containing input fastq files.
#' @param dir_out Output directory. If it does not exist, it will be created.
#' @param primer_16S_fwd,primer_16S_rev DNA sequences of 16S forward and reverse primer, respectively
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#' @param post_samplename_pattern1,post_samplename_pattern2 (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#'
#' @return Integer matrix denoting the number of reads remaining after primer-trimming for each input file.
#' @export
#'
#' @examples
#' \dontrun{
#' trimPrimers16S(c("sample1_R1.fastq", "sample1_R2.fastq", "sample2_R1.fastq", "sample2_R2.fastq"), "path/to/input", "path/to/output", "CCTACGGGNBGCASCAG", "GACTACNVGGGTATCTAATCC", multithread = TRUE)
#' }
trimPrimers16S <- function(fn, dir_in, dir_out, primer_16S_fwd, primer_16S_rev, multithread = MULTITHREAD, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq"#, quiet=T
){
  fn_fullname <- file.path(dir_in, fn)

  fnFs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern1, fn_fullname)]
  fnRs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern2, fn_fullname)]
  if(length(fnFs) + length(fnRs) == 0) warning(paste0("trimPrimer16S: ", "No files found at specified location(s) within ", dir_in, ". Check file path, or post_samplename_pattern argument(s)."))

  # Extract sample names
  sample.names <- getSampleName(fnFs)

  # Create output directory
  if(!dir.exists(dir_out)) dir.create(dir_out)

  # Create primer-trimmed file paths in output directory
  fnFs.cut <- file.path(dir_out, basename(fnFs))
  fnRs.cut <- file.path(dir_out, basename(fnRs))

  # Just trim the length of fwd and rev primers
  qa.out <- filterAndTrim(fnFs, fnFs.cut, fnRs, fnRs.cut, multithread = multithread,
                          matchIDs = TRUE, trimLeft = c(nchar(primer_16S_fwd), nchar(primer_16S_rev)), compress=F)

  colnames(qa.out) <- c("input", "trimmed")
  rownames(qa.out) <- sample.names
  return(qa.out)
}


#' Trim Primers from 16S Sequences (with metadata)
#'
#' Trims primers from 16S sequences using \code{\link[dada2]{filterAndTrim}}. This function
#' assumes that each read begins with its full primer sequence and operates
#' by truncating the beginning of each read by the length of its primer
#' sequence.
#'
#' @param fn Full names of input fastq files, including directory. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param dir_out Directory where filtered fastq files will be written.
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param primer_16S_fwd,primer_16S_rev DNA sequences of 16S forward and reverse primer, respectively
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#'
#' @return Integer matrix denoting the number of reads remaining after primer-trimming for each input file.
#' @export
#'
#' @examples
#' \dontrun{
#' trimPrimers16S(c("sample1_R1.fastq", "sample1_R2.fastq", "sample2_R1.fastq", "sample2_R2.fastq"), "path/to/output", meta, "CCTACGGGNBGCASCAG", "GACTACNVGGGTATCTAATCC", multithread = TRUE)
#' }
trimPrimers16S2 <- function(fn, dir_out, meta, primer_16S_fwd, primer_16S_rev, multithread = MULTITHREAD) {
  dir.create(dir_out, recursive = TRUE)
  fn_out <- file.path(dir_out, basename(fn))

  keep_fn <- file.exists(fn) & !duplicated(fn)
  fn <- fn[keep_fn]
  if(length(fn) == 0) warning(paste0("trimPrimers16S: No files found at specified location(s). Check file paths, or input metadata."))

  # Get metadata matching files
  meta_ext <- matchFastqToMetadata(fn, meta)

  # Reference metadata to retrieve R1 and R2 files
  fn_pairs <- getPairedFastqFiles(fn, meta_ext[,-1], value=FALSE)
  fnFs <- fn[fn_pairs[[1]]]
  fnRs <- fn[fn_pairs[[2]]]

  # Confirm target gene
  if(any(!grepl("16S", meta_ext$targetGene))) warning("You are using trimPrimers16S() on some non-16S files. Did you mean to use trimPrimersITS()?")

  # Attempt to get DNA sample IDs to use as row names
  dnaSampleIDs <- as.character(meta_ext$dnaSampleID[keep_fn][fn_pairs[[1]]])
  fn_as_rownames <- any(is.na(dnaSampleIDs))
  if(fn_as_rownames) {
    warning(sum(is.na(dnaSampleIDs)),
            " files do not have associated dnaSampleIDs in the input metadata. ",
            "Therefore, input R1 file names, instead of dnaSampleIDs, will be used as the rownames for the read tracking table.")
  }

  # Confirm output filenames
  trimFs <- fn_out[keep_fn][fn_pairs[[1]]]
  trimRs <- fn_out[keep_fn][fn_pairs[[2]]]

  # Just trim the length of fwd and rev primers
  qa.out <- filterAndTrim(fnFs, trimFs, fnRs, trimRs, multithread = multithread,
                          matchIDs = TRUE, trimLeft = c(nchar(primer_16S_fwd), nchar(primer_16S_rev)), compress=F)

  colnames(qa.out) <- c("input", "trimmed")
  if(fn_as_rownames) {
    rownames(qa.out) <- fnFs
  } else {
    rownames(qa.out) <- dnaSampleIDs
  }
  return(qa.out)
}



#' Trim Primers from ITS Sequences
#'
#' Trims primers from ITS sequences using cutadapt. Cutadapt must be installed in order for this to work. Currently only supports R1 (forward-read) files.
#'
#' @param fn Names of input fastq files, excluding directory path which is specified by dir_in. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning. It is assumed that these are R1 (forward-read) files only.
#' @param dir_in Directory containing input fastq files.
#' @param dir_out Output directory. If it does not exist, it will be created.
#' @param primer_ITS_fwd,primer_ITS_rev DNA sequence of the ITS forward and reverse primer, respectively.
#' @param cutadapt_path Default CUTADAPT_PATH in params.R. Path to cutadapt on your file system.
#' @param post_samplename_pattern1,post_samplename_pattern2 (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#' @param very_verbose Default FALSE. Whether to print output from cutadapt. Unlike some other "verbose" arguments associated with the functions in this package, this does not default to VERBOSE in params.R.
#' @param discard_untrimmed Default FALSE. Whether to discard reads where a primer could not be found, leaving only those reads in which a primer has been trimmed.
#'
#' @return No value is returned.
#' @export
#'
#' @examples
#' \dontrun{
#' trimPrimersITS(c("sample1_ITS_R1.fastq", "sample2_ITS_R1.fastq"), "path/to/input", "path/to/output", "CTTGGTCATTTAGAGGAAGTAA")
#' }
trimPrimersITS <- function(fn, dir_in, dir_out, primer_ITS_fwd, primer_ITS_rev, cutadapt_path = CUTADAPT_PATH, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq", very_verbose=FALSE, discard_untrimmed=FALSE#, quiet=T
){

  fn_fullname <- file.path(dir_in, fn)

  fnFs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern1, fn_fullname)]
  if(length(fnFs) == 0) warning(paste0("trimPrimerITS: ", "No files found at specified location(s) within ", dir_in, ". Check file path, or post_samplename_pattern argument(s)."))

  # Extract sample names
  sample.names <- getSampleName(fnFs)

  # Create output directory
  if(!dir.exists(dir_out)) dir.create(dir_out)

  # Create primer-trimmed file paths in output directory
  fnFs.cut <- file.path(dir_out, basename(fnFs))

  # Get reverse-complement sequence of reverse primer
  rev_rc <- dada2:::rc(primer_ITS_rev)

  # Trim forward primer and the reverse-complement of reverse primer off of R1 (forward reads)
  R1.flags <- paste("-g", primer_ITS_fwd, "-a", rev_rc)

  # Discard untrimmed reads?
  if(discard_untrimmed) {
    discard_untrimmed_flag <- "--discard-untrimmed"
  } else {
    discard_untrimmed_flag <- ""
  }

  # Run Cutadapt
  for(i in seq_along(fnFs)) {
    system2(cutadapt_path, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                                    "-o", fnFs.cut[i], # output files
                                    fnFs[i], # input files
                                    "--minimum-length", "1", # min length of cutadapted reads: >0
                                    discard_untrimmed_flag,
                                    "-e", "0.2"), # -e 0.2 allows up to 4 mismatched bases
            stdout = ifelse(very_verbose, "", FALSE))
  }
}


#' Trim Primers from ITS Sequences (with metadata)
#'
#' Trims primers from ITS sequences using cutadapt. Cutadapt must be installed in order for this to work. Currently only supports forward-read (R1) sequences. If reverse reads are included, they will be ignored.
#'
#' @param fn Full names of input fastq files, including directory. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning. It is assumed that these are R1 (forward-read) files only.
#' @param dir_out Directory where filtered fastq files will be written.
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param primer_ITS_fwd,primer_ITS_rev DNA sequence of the ITS forward and reverse primer, respectively.
#' @param cutadapt_path Default CUTADAPT_PATH in params.R. Path to cutadapt on your file system.
#' @param very_verbose Default FALSE. Whether to print output from cutadapt. Unlike some other "verbose" arguments associated with the functions in this package, this does not default to VERBOSE in params.R.
#' @param discard_untrimmed Default FALSE. Whether to discard reads where a primer could not be found, leaving only those reads in which a primer has been trimmed.
#'
#' @return No value is returned.
#' @export
#'
#' @examples
#' \dontrun{
#' trimPrimersITS(c("sample1_ITS_R1.fastq", "sample2_ITS_R1.fastq"), "path/to/input", "path/to/output", "CTTGGTCATTTAGAGGAAGTAA")
#' }
trimPrimersITS2 <- function(fn, dir_out, meta, primer_ITS_fwd, primer_ITS_rev, cutadapt_path = CUTADAPT_PATH, very_verbose=FALSE, discard_untrimmed=FALSE) {
  dir.create(dir_out, recursive = TRUE)
  fn_out <- file.path(dir_out, basename(fn))

  keep_fn <- file.exists(fn) & !duplicated(fn)
  fn <- fn[keep_fn]
  if(length(fn) == 0) warning(paste0("trimPrimersITS: No files found at specified location(s). Check file paths, or input metadata."))

  # Reference metadata to retrieve R1 files
  meta_ext <- matchFastqToMetadata(fn, meta)
  fnFs <- fn[grep("R1", meta_ext$rawDataFileDescription)]

  # Confirm target gene
  if(any(!grepl("ITS", meta_ext$targetGene))) warning("You are using trimPrimersITS() on some non-ITS files. Did you mean to use trimPrimers16S()?")

  # Attempt to get DNA sample IDs to use as row names
  dnaSampleIDs <- as.character(meta_ext$dnaSampleID[match(fnFs, meta_ext$file)])
  fn_as_rownames <- any(is.na(dnaSampleIDs))
  if(fn_as_rownames) {
    warning(sum(is.na(dnaSampleIDs)),
            " files do not have associated dnaSampleIDs in the input metadata. ",
            "Therefore, input R1 file names, instead of dnaSampleID, will be used as the rownames for the read tracking table.")
  }

  # Confirm output filenames
  trimFs <- fn_out[keep_fn][grep("R1", meta_ext$rawDataFileDescription)]

  # Get reverse-complement sequence of reverse primer
  rev_rc <- dada2:::rc(primer_ITS_rev)

  # Trim forward primer and the reverse-complement of reverse primer off of R1 (forward reads)
  R1.flags <- paste("-g", primer_ITS_fwd, "-a", rev_rc)

  # Discard untrimmed reads?
  if(discard_untrimmed) {
    discard_untrimmed_flag <- "--discard-untrimmed"
  } else {
    discard_untrimmed_flag <- ""
  }

  # Run Cutadapt
  for(i in seq_along(fnFs)) {
    system2(cutadapt_path, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                                    "-o", trimFs[i], # output files
                                    fnFs[i], # input files
                                    "--minimum-length", "1", # min length of cutadapted reads: >0
                                    discard_untrimmed_flag,
                                    "-e", "0.2"), # -e 0.2 allows up to 4 mismatched bases
            stdout = ifelse(very_verbose, "", FALSE))
  }
}


#' Filter 16S Sequences
#'
#' Applies a quality filter to 16S sequence fastq files via the \code{\link[dada2]{filterAndTrim}} function.
#' It is assumed that both forward- and reverse-read files are included.
#'
#' @param fn Names of input fastq files, excluding directory path which is specified by dir_in. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param dir_in Directory containing input fastq files.
#' @param dir_out Path to output directory where filtered fastq files will be written.
#' @param trunc_qscore Default 23. Quality score at which point to truncate each read, if truncLen is NULL.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#' @param post_samplename_pattern1,post_samplename_pattern2 (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#' @param ... Other arguments to be passed to \code{\link[dada2]{filterAndTrim}}, such as maxEE. See documentation for more details.
#'
#' @return Two-column matrix displaying the number of reads in input vs. output for each file.
#' @export
qualityFilter16S <- function(fn, dir_in, dir_out, trunc_qscore = 23, multithread = MULTITHREAD, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq", ...){
  fn_fullname <- file.path(dir_in, fn)

  fnFs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern1, fn_fullname)]
  fnRs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern2, fn_fullname)]
  if(length(fnFs) + length(fnRs) == 0) warning(paste0("qualityFilter16S: ", "No files found at specified location(s) within ", dir_in, ". Check file path, or post_samplename_pattern argument(s)."))

  sample.names <- getSampleName(fnFs) # extract sample names
  filtFs <- file.path(dir_out, basename(fnFs)) # create filtered filenames
  filtRs <- file.path(dir_out, basename(fnRs)) # create filtered filenames

  dots <- list(...)

  # If truncLen is NOT in user-provided arguments
  if (!("truncLen" %in% names(dots))) {
    # GETTING TRUNCATION LENGTH USING A SUBSET OF QUALITY SCORES
    n_files <- min(length(unique(sample.names)), 30)

    fwd.trunc.lengths <- list()
    rev.trunc.lengths <- list()
    fwd.trunc.lengths[1:n_files] <- getTruncationLength(fnFs[1:n_files], verbose = VERBOSE,
                                                        qscore = trunc_qscore)
    rev.trunc.lengths[1:n_files] <- getTruncationLength(fnRs[1:n_files], verbose = VERBOSE,
                                                        qscore = trunc_qscore)

    # remove any samples that have low-quality reads early on, to avoid spoiling the whole run
    # idk how to loop this command
    to_remove <- which(unlist(fwd.trunc.lengths) < 100 | unlist(rev.trunc.lengths) < 100)
    if (length(to_remove) > 0){
      fwd.trunc.lengths <- fwd.trunc.lengths[-to_remove]
      rev.trunc.lengths <- rev.trunc.lengths[-to_remove]
    }

    # if (mean == T){
    #   fwd.trunc.length <- round(mean(unlist(fwd.trunc.lengths))) # tested this with mean vs minimum quality-score threshold;
    #   rev.trunc.length <- round(mean(unlist(rev.trunc.lengths))) # minimum retained more reads after filtering.
    # } else {
    fwd.trunc.length <- round(min(unlist(fwd.trunc.lengths))) # tested this with mean vs minimum quality-score threshold;
    rev.trunc.length <- round(min(unlist(rev.trunc.lengths))) # minimum retained more reads after filtering.
    # }

    # Set the minimum lengths for fwd/reverse reads (if they're too short, they cannot be merged).
    if (rev.trunc.length < 200) rev.trunc.length <- 200
    if (fwd.trunc.length < 245) fwd.trunc.length <- 245
    cat(paste0("Fwd truncation length: ", fwd.trunc.length, "\nRev truncation length: ", rev.trunc.length, "\n"))
    truncLen <- c(fwd.trunc.length, rev.trunc.length)

    # Else, if truncLen is in user-provided arguments
  } else {
    message("Using truncLen ", paste(dots$truncLen, collapse=" "))
    arguments <- c(list(fnFs, filtFs, fnRs, filtRs, compress=TRUE, multithread=multithread), dots)
  }

  out <- do.call(filterAndTrim, arguments)

  rownames(out) <- sample.names

  return(out)
}
# Removed this: #' @param truncLen_manual Default NULL. Single integer: truncation length to use across all files. Two-integer vector: truncation length to use for the forward-read and reverse-read files, respectively. If NULL (default), determines truncation length(s) based on \code{\link{getTruncationLength}} with a quality score threshold of trunc_qscore.


#' Filter 16S Sequences (with metadata)
#'
#' Applies a quality filter to 16S sequence fastq files via the \code{\link[dada2]{filterAndTrim}} function.
#' It is assumed that both forward- and reverse-read files are included.
#'
#' @param fn Full names of input fastq files, including directory. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning. It is assumed that these are R1 (forward-read) files only.
#' @param dir_out Directory where filtered fastq files will be written.
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param trunc_qscore Default 23. Quality score at which point to truncate each read, if truncLen is NULL.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#' @param ... Other arguments to be passed to \code{\link[dada2]{filterAndTrim}}, such as maxEE and truncLen. See documentation for more details.
#'
#' @return Two-column matrix displaying the number of reads in input vs. output for each file.
#' @export
qualityFilter16S2 <- function(fn, dir_out, meta, trunc_qscore = 23, multithread = MULTITHREAD, ...){
  dir.create(dir_out, recursive = TRUE)
  fn_out <- file.path(dir_out, basename(fn))

  keep_fn <- file.exists(fn) & !duplicated(fn)
  fn <- fn[keep_fn]
  if(length(fn) == 0) warning(paste0("qualityFilter16S: No files found at specified location(s). Check file paths, or input metadata."))

  # Get metadata matching files
  meta_ext <- matchFastqToMetadata(fn, meta)

  # Reference metadata to retrieve R1 and R2 files
  fn_pairs <- getPairedFastqFiles(fn, meta_ext[,-1], value=FALSE)
  fnFs <- fn[fn_pairs[[1]]]
  fnRs <- fn[fn_pairs[[2]]]

  # Confirm target gene
  if(any(!grepl("16S", meta_ext$targetGene))) warning("You are using qualityFilter16S() on some non-16S files. Did you mean to use qualityFilterITS()?")

  # Attempt to get DNA sample IDs to use as row names
  dnaSampleIDs <- as.character(meta_ext$dnaSampleID[keep_fn][fn_pairs[[1]]])
  fn_as_rownames <- any(is.na(dnaSampleIDs))
  if(fn_as_rownames) {
    warning(sum(is.na(dnaSampleIDs)),
            " files do not have associated dnaSampleIDs in the input metadata. ",
            "Therefore, input R1 file names, instead of dnaSampleIDs, will be used as the rownames for the read tracking table.")
  }

  # Confirm output filenames
  filtFs <- fn_out[keep_fn][fn_pairs[[1]]]
  filtRs <- fn_out[keep_fn][fn_pairs[[2]]]

  # Read in dots args
  dots <- list(...)

  # If truncLen is NOT in user-provided arguments
  if (!("truncLen" %in% names(dots))) {
    # GETTING TRUNCATION LENGTH USING A SUBSET OF QUALITY SCORES
    n_files <- min(length(fnFs), 30)

    fwd.trunc.lengths <- list()
    rev.trunc.lengths <- list()
    fwd.trunc.lengths[1:n_files] <- getTruncationLength(fnFs[1:n_files], verbose = VERBOSE,
                                                        qscore = trunc_qscore)
    rev.trunc.lengths[1:n_files] <- getTruncationLength(fnRs[1:n_files], verbose = VERBOSE,
                                                        qscore = trunc_qscore)

    # remove any samples that have low-quality reads early on, to avoid spoiling the whole run
    # idk how to loop this command
    to_remove <- which(unlist(fwd.trunc.lengths) < 100 | unlist(rev.trunc.lengths) < 100)
    if (length(to_remove) > 0){
      fwd.trunc.lengths <- fwd.trunc.lengths[-to_remove]
      rev.trunc.lengths <- rev.trunc.lengths[-to_remove]
    }

    # if (mean == T){
    #   fwd.trunc.length <- round(mean(unlist(fwd.trunc.lengths))) # tested this with mean vs minimum quality-score threshold;
    #   rev.trunc.length <- round(mean(unlist(rev.trunc.lengths))) # minimum retained more reads after filtering.
    # } else {
    fwd.trunc.length <- round(min(unlist(fwd.trunc.lengths))) # tested this with mean vs minimum quality-score threshold;
    rev.trunc.length <- round(min(unlist(rev.trunc.lengths))) # minimum retained more reads after filtering.
    # }

    # Set the minimum lengths for fwd/reverse reads (if they're too short, they cannot be merged).
    if (rev.trunc.length < 200) rev.trunc.length <- 200
    if (fwd.trunc.length < 245) fwd.trunc.length <- 245
    message("Using auto-selected truncLen:")
    cat(paste0("Fwd truncation length: ", fwd.trunc.length, "\nRev truncation length: ", rev.trunc.length, "\n"))
    dots$truncLen <- c(fwd.trunc.length, rev.trunc.length)

    # Else, if truncLen is in user-provided arguments
  } else {
    message("Using user-provided truncLen: ", paste(dots$truncLen, collapse=" "))
  }

  # Run dada::filterAndTrim
  arguments <- c(list(fnFs, filtFs, fnRs, filtRs, compress=TRUE, multithread=multithread), dots)
  out <- do.call(filterAndTrim, arguments)

  if(fn_as_rownames) {
    rownames(out) <- fnFs
  } else {
    rownames(out) <- dnaSampleIDs
  }

  return(out)
}


#' Filter ITS Sequences
#'
#' Applies a quality filter to ITS sequence fastq files via the \code{\link[dada2]{filterAndTrim}} function.
#' Currently only supports filtering forward-read sequences.
#'
#' @param fn Names of input fastq files, excluding directory path which is specified by dir_in. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning. It is assumed that these are R1 (forward-read) files only.
#' @param dir_in Directory containing input fastq files.
#' @param dir_out Path to output directory where filtered fastq files will be written.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#' @param post_samplename_pattern1,post_samplename_pattern2 (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#' @param ... Other arguments to be passed to \code{\link[dada2]{filterAndTrim}}, such as maxEE. See documentation for more details.
#'
#' @return Two-column matrix displaying the number of reads in input vs. output for each file.
#' @export
qualityFilterITS <- function(fn, dir_in, dir_out, multithread = MULTITHREAD, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq", ...){

  fn_fullname <- file.path(dir_in, fn)

  fnFs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern1, fn_fullname)]
  if(length(fnFs) == 0) warning(paste0("trimPrimerITS: ", "No files found at specified location(s) within ", dir_in, ". Check file path, or post_samplename_pattern argument(s)."))

  sample.names <- getSampleName(fnFs) # extract sample names

  gz_ind <- grepl("\\.gz$", basename(fnFs))
  if(!all(gz_ind)) {
    filtFs <- file.path(dir_out, basename(fnFs))
    filtFs[!gz_ind] <- paste0(filtFs, ".gz")
  } else {
    filtFs <- file.path(dir_out, basename(fnFs))
  }

  out <- filterAndTrim(fnFs, filtFs,
                       compress = TRUE,
                       matchIDs = TRUE,
                       multithread = multithread,
                       ...)

  rownames(out) <- sample.names

  return(out)
}

#' Filter ITS Sequences (with metadata)
#'
#' Applies a quality filter to ITS sequence fastq files via the \code{\link[dada2]{filterAndTrim}} function.
#' Currently only supports filtering forward-read (R1) sequences. If reverse reads are included, they will be ignored.
#'
#' @param fn Full names of input fastq files, including directory. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning. It is assumed that these are R1 (forward-read) files only.
#' @param dir_out Directory where filtered fastq files will be written.
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#' @param ... Other arguments to be passed to \code{\link[dada2]{filterAndTrim}}, such as maxEE. See documentation for more details.
#'
#' @return Two-column matrix displaying the number of reads in input vs. output for each file.
#' @export
qualityFilterITS2 <- function(fn, dir_out, meta, multithread = MULTITHREAD, ...){
  dir.create(dir_out, recursive = TRUE)
  fn_out <- file.path(dir_out, basename(fn))

  keep_fn <- file.exists(fn) & !duplicated(fn)
  fn <- fn[keep_fn]
  if(length(fn) == 0) warning(paste0("qualityFilterITS: No files found at specified location(s). Check file paths, or input metadata."))

  # Reference metadata to retrieve R1 files
  meta_ext <- matchFastqToMetadata(fn, meta)
  fnFs <- fn[grep("R1", meta_ext$rawDataFileDescription)]

  # Confirm target gene
  if(any(!grepl("ITS", meta_ext$targetGene))) warning("You are using qualityFilterITS() on some non-ITS files. Did you mean to use qualityFilter16S()?")

  # Attempt to get DNA sample IDs to use as row names
  dnaSampleIDs <- as.character(meta_ext$dnaSampleID[match(fnFs, meta_ext$file)])
  fn_as_rownames <- any(is.na(dnaSampleIDs))
  if(fn_as_rownames) {
    warning(sum(is.na(dnaSampleIDs)),
            " files do not have associated dnaSampleIDs in the input metadata. ",
            "Therefore, input R1 file names, instead of dnaSampleID, will be used as the rownames for the read tracking table.")
  }

  # Confirm output filenames
  filtFs <- fn_out[keep_fn][grep("R1", meta_ext$rawDataFileDescription)]

  # Run quality filter
  out <- filterAndTrim(fnFs, filtFs,
                       compress = TRUE,
                       multithread = multithread,
                       ...)

  if(fn_as_rownames) {
    rownames(out) <- fnFs
  } else {
    rownames(out) <- dnaSampleIDs
  }
  return(out)
}

#' Run Dada on Paired-End 16S Sequences
#'
#' Runs the Dada algorithm to infer sample composition from paired-end 16S fastq files.
#' This implementation is based on Ben Callahan's vignette at \url{https://benjjneb.github.io/dada2/bigdata_paired.html}.
#' @param fn Names of input fastq files, excluding directory path which is specified by dir_in. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param dir_in Directory containing input fastq files.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading.
#' @param VERBOSE Default FALSE. Whether to print messages regarding the dimensions of the resulting sequence table and the distribution of sequence lengths.
#' @param seed (Optional) Integer to use as random seed for reproducibility.
#' @param nbases (Optional) Number of bases to use for learning errors. Default 1e7.
#'
#' @return A list of three elements. \strong{seqtab} is the sequence table before removing chimeras, \strong{seqtab.nochim} is the sequence table after removing chimeras, and \strong{track} is a data frame displaying the number of reads remaining for each sample at various points throughout the processing pipeline.
#' @export
#'
#' @examples
#' \dontrun{
#' seqtab.list <- runDada16S(c("sample1_R1.fastq", "sample1_R2.fastq", "sample2_R1.fastq", "sample2_R2.fastq"), seed=1010100)
#' }
runDada16S <- function(fn, dir_in, multithread = MULTITHREAD, verbose = FALSE, seed = NULL, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq", nbases=1e7){
  if (!is.null(seed)) set.seed(seed)

  fn_fullname <- file.path(dir_in, fn)

  filtFs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern1, fn_fullname)]
  filtRs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern2, fn_fullname)]
  if(length(filtFs) + length(filtRs) == 0) warning(paste0("runDada16S: ", "No files found at specified location(s) within ", dir_in, ". Check file path, or post_samplename_pattern argument(s)."))

  # Keep files with counterpart
  matched_files <- remove_unmatched_files(filtFs, filtRs)
  filtFs <- matched_files[["R1"]]
  filtRs <- matched_files[["R2"]]

  # Create outnames
  sample.names <- sapply(strsplit(basename(filtFs), paste0("(",post_samplename_pattern1,")|(",post_samplename_pattern2,")")), `[`, 1)
  sample.namesR <- sapply(strsplit(basename(filtRs), paste0("(",post_samplename_pattern1,")|(",post_samplename_pattern2,")")), `[`, 1)
  # sample.names <- sapply(strsplit(basename(filtFs), "_R(1|2).fastq"), `[`, 1) # Assumes filename = samplename_RX.fastq.gz
  # sample.namesR <- sapply(strsplit(basename(filtRs), "_R(1|2).fastq"), `[`, 1) # Assumes filename = samplename_RX.fastq.gz
  if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  # Learn forward and reverse error rates
  errF <- learnErrors(filtFs, nbases=nbases, randomize=TRUE, multithread=multithread, verbose=verbose)
  errR <- learnErrors(filtRs, nbases=nbases, randomize=TRUE, multithread=multithread, verbose=verbose)

  # Create output vectors
  derepF <- derepR <- ddF <- ddR <- mergers <- vector("list", length(sample.names))
  names(derepF) <- names(derepR) <- names(ddF) <- names(ddR) <- names(mergers) <- sample.names

  # Sample inference and merger of paired-end reads
  for(i in 1:length(sample.names)) {
    sam <- sample.names[[i]]
    cat("Processing:", sam, "\n")

    derepF[[i]] <- derepFastq(filtFs[[i]])
    derepR[[i]] <- derepFastq(filtRs[[i]])

    ddF[[i]] <- dada(derepF[[i]], err=errF, multithread=multithread, verbose=verbose)
    ddR[[i]] <- dada(derepR[[i]], err=errR, multithread=multithread, verbose=verbose)
    mergers[[i]] <- mergePairs(ddF[[i]], derepF[[i]], ddR[[i]], derepR[[i]], maxMismatch=1, minOverlap = 6, verbose=verbose)
    cat(paste("Exact sequence variants inferred for sample:", sam,". \n"))
  }

  # Construct sequence table and remove chimeras
  seqtab <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=multithread, verbose=verbose)

  if(verbose){
    cat(paste0("\n\nDimensions of ESV table: ", dim(seqtab)[1], " samples, ", dim(seqtab)[2], " ESVs\n"))
    cat("\nDistribution of sequence lengths (should be bimodal for v3v4 region):")
    print(table(nchar(getSequences(seqtab))))
    cat(paste0("\n", round(sum(seqtab.nochim)/sum(seqtab), 3)*100, "% reads remain after removal of chimeras"))
  }

  # Track reads through this part of the pipeline
  suppressWarnings({
    track <- Reduce(
      function(x, y, ...) transform(merge(x, y, by=0, all = TRUE, ...), row.names=Row.names, Row.names=NULL),
      list(sapply(derepF, getN),
           sapply(derepR, getN),
           sapply(ddF, getN),
           denoisedR = sapply(ddR, getN),
           sapply(mergers, getN),
           rowSums(seqtab.nochim))
    )
  })
  names(track) <- c("derepF", "derepR", "denoisedF", "denoisedR", "merged", "nonchim")
  track[is.na(track)] <- 0

  return(list("seqtab" = seqtab,
              "seqtab.nochim" = seqtab.nochim,
              "track" = track))
}


#' Run Dada on Paired-End 16S Sequences (with metadata)
#'
#' Runs the Dada algorithm to infer sample composition from paired-end 16S fastq files.
#' This implementation is based on Ben Callahan's vignette at \url{https://benjjneb.github.io/dada2/bigdata_paired.html}.
#'
#' @param fn Full names of input fastq files, including directory. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param out_seqtab,out_track (Optional) File locations where sequence table and read-tracking table (respectively) will be saved as csv files. If blank (default), will not save to file.
#' @param remove_chimeras (Optional) Default TRUE. Whether to remove chimeras from the sequence table. Currently only supports removal using \code{\link[dada2]{removeBimeraDenovo}} with the consensus method.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading.
#' @param verbose Default FALSE. Whether to print messages regarding the samples being processed, dimensions of the resulting sequence table, and the distribution of sequence lengths.
#' @param seed (Optional) Integer to use as random seed for reproducibility.
#' @param nbases (Optional) Number of bases to use for learning errors. Default 1e7.
#'
#' @return A list of two elements. \strong{seqtab} is the sequence table, with chimeras removed if remove_chimeras == TRUE, and \strong{track} is a data frame displaying the number of reads remaining for each sample at various points throughout the processing pipeline.
#' @export
#'
#' @examples
#' \dontrun{
#' seqtab.list <- runDada16S(c("sample1_R1.fastq", "sample1_R2.fastq", "sample2_R1.fastq", "sample2_R2.fastq"), meta, seed=1010100)
#' }
runDada16S2 <- function(fn, meta, out_seqtab = "", out_track = "", remove_chimeras = TRUE, multithread = MULTITHREAD, verbose = FALSE, seed = NULL, nbases=1e7){
  if (!is.null(seed)) set.seed(seed)

  keep_fn <- file.exists(fn) & !duplicated(fn)
  fn <- fn[keep_fn]
  if(length(fn) == 0) warning(paste0("runDada16S: No files found at specified location(s). Check file paths, or input metadata."))

  # Get metadata matching files
  meta_ext <- matchFastqToMetadata(fn, meta)

  # Reference metadata to retrieve R1 and R2 files
  fn_pairs <- getPairedFastqFiles(fn, meta_ext[,-1], value=FALSE)
  fnFs <- fn[fn_pairs[[1]]]
  fnRs <- fn[fn_pairs[[2]]]

  # Confirm target gene
  if(any(!grepl("16S", meta_ext$targetGene))) warning("You are using runDada16S() on some non-16S files. Did you mean to use runDadaITS()?")

  # Attempt to get DNA sample IDs to use as row names
  dnaSampleIDs <- as.character(meta_ext$dnaSampleID[match(fnFs, meta_ext$file)])
  fn_as_rownames <- any(is.na(dnaSampleIDs))
  if(fn_as_rownames) {
    warning(sum(is.na(dnaSampleIDs)),
            " files do not have associated dnaSampleIDs in the input metadata. ",
            "Therefore, input R1 file names, instead of dnaSampleIDs, will be used as the rownames for the read tracking table.")
  }

  names(fnFs) <- dnaSampleIDs
  names(fnRs) <- dnaSampleIDs

  # Learn forward and reverse error rates
  errF <- learnErrors(fnFs, nbases=nbases, randomize=TRUE, multithread=multithread, verbose=verbose)
  errR <- learnErrors(fnRs, nbases=nbases, randomize=TRUE, multithread=multithread, verbose=verbose)

  # Create output vectors
  derepF <- derepR <- ddF <- ddR <- mergers <- vector("list", length(dnaSampleIDs))
  names(derepF) <- names(derepR) <- names(ddF) <- names(ddR) <- names(mergers) <- dnaSampleIDs

  # Sample inference and merger of paired-end reads
  message("runDada16S progress out of ", length(dnaSampleIDs), " samples:")
  progressbar <- txtProgressBar(min = 0, max = length(dnaSampleIDs), style = 3)
  for(i in 1:length(dnaSampleIDs)) {
    sam <- dnaSampleIDs[[i]]
    if(verbose) cat("Processing:", sam, "\n")

    derepF[[i]] <- derepFastq(fnFs[[i]])
    derepR[[i]] <- derepFastq(fnRs[[i]])

    ddF[[i]] <- dada(derepF[[i]], err=errF, multithread=multithread, verbose=verbose)
    ddR[[i]] <- dada(derepR[[i]], err=errR, multithread=multithread, verbose=verbose)
    if(verbose) {
      mergers[[i]] <- mergePairs(ddF[[i]], derepF[[i]], ddR[[i]], derepR[[i]], maxMismatch=1, minOverlap = 6, verbose=TRUE)
    } else {
      mergers[[i]] <- suppressMessages(mergePairs(ddF[[i]], derepF[[i]], ddR[[i]], derepR[[i]], maxMismatch=1, minOverlap = 6, verbose=FALSE))
    }
    if(verbose) cat(paste("Exact sequence variants inferred for sample:", sam,". \n"))
    setTxtProgressBar(progressbar, i)
  }
  close(progressbar)

  # Construct sequence table and remove chimeras
  if(verbose) {
    seqtab <- makeSequenceTable(mergers)
  } else {
    seqtab <- suppressMessages(makeSequenceTable(mergers))
  }
  if(remove_chimeras) seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=multithread, verbose=verbose)

  if(verbose){
    cat(paste0("\nDimensions of ESV table: ", dim(seqtab)[1], " samples, ", dim(seqtab)[2], " ESVs\n"))
    cat("\nDistribution of sequence lengths (should be bimodal for v3v4 region):")
    print(table(nchar(getSequences(seqtab))))
    if(remove_chimeras) cat(paste0("\n", round(sum(seqtab.nochim)/sum(seqtab), 3)*100, "% reads remain after removal of chimeras"))
  }

  # Track reads through this part of the pipeline
  suppressWarnings({
    track <- Reduce(
      function(x, y, ...) transform(merge(x, y, by=0, all = TRUE, ...), row.names=Row.names, Row.names=NULL),
      if(remove_chimeras) {
        list(sapply(derepF, getN),
             sapply(derepR, getN),
             sapply(ddF, getN),
             denoisedR = sapply(ddR, getN),
             suppressMessages(sapply(mergers, getN)),
             rowSums(seqtab.nochim))
      } else {
        list(sapply(derepF, getN),
             sapply(derepR, getN),
             sapply(ddF, getN),
             denoisedR = sapply(ddR, getN),
             suppressMessages(sapply(mergers, getN)))
      }
    )
  })
  if(remove_chimeras) {
    names(track) <- c("derepF", "derepR", "denoisedF", "denoisedR", "merged", "nonchim")
  } else {
    names(track) <- c("derepF", "derepR", "denoisedF", "denoisedR", "merged")
  }
  track[is.na(track)] <- 0

  if(out_seqtab != "") {
    if(!dir.exists(dirname(out_seqtab))) {
      dir.create(dirname(out_seqtab), recursive=TRUE)
      message("Output directory created for sequence table: ", dirname(out_seqtab))
    }
    if(remove_chimeras) {
      write.csv(seqtab.nochim, out_seqtab)
    } else {
      write.csv(seqtab, out_seqtab)
    }
  }
  if(out_track != "") {
    if(!dir.exists(dirname(out_track))) {
      dir.create(dirname(out_track), recursive=TRUE)
      message("Output directory created for read tracking table: ", dirname(out_track))
    }
    write.csv(track, out_track)
  }

  if(remove_chimeras) {
    return(list("seqtab" = seqtab.nochim,
                "track" = track))
  } else {
    return(list("seqtab" = seqtab,
                "track" = track))
  }
}


#' Run Dada on R1 ITS sequences
#'
#' Runs the Dada algorithm to infer sample composition from R1 ITS fastq files.
#' We recommend using only the R1 reads because attempting to merge with R2
#' reads may result in less accurate representations of the fungal community
#' composition (Pauvert et al., 2019).
#' This implementation is based on Ben Callahan's vignettes at \url{https://benjjneb.github.io/dada2/bigdata.html}
#' and \url{https://benjjneb.github.io/dada2/ITS_workflow.html}.
#'
#' @param fn Names of input fastq files, excluding directory path which is specified by dir_in. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning. It is assumed that these are R1 (forward-read) files only.
#' @param dir_in Directory containing input fastq files.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading.
#' @param verbose Default FALSE. Whether to print messages regarding the dereplication step, the denoising step, and the dimensions of the resulting sequence table and the distribution of sequence lengths.
#' @param seed (Optional) Integer to use as random seed for reproducibility.
#' @param post_samplename_pattern1,post_samplename_pattern2 (Optional) Character pattern within the filename which immediately follows the end of the sample name. Defaults to "_R(1|2).*\\.fastq", as NEON fastq files typically consist of a sample name followed by "_R1.fastq" or "_R2.fastq", etc.
#' @param nbases (Optional) Number of bases to use for learning errors. Default 1e7.
#'
#' @return A list of three elements. \strong{seqtab} is the sequence table before removing chimeras, \strong{seqtab.nochim} is the sequence table after removing chimeras, and \strong{track} is a data frame displaying the number of reads remaining for each sample at various points throughout the processing pipeline.
#' @export
#'
#' @examples
#' \dontrun{
#' seqtab.list <- runDadaITS(c("sample1_R1.fastq", "sample1_R2.fastq", "sample2_R1.fastq", "sample2_R2.fastq"), './seq/filtered/', seed=1010100)
#' }
runDadaITS <- function(fn, dir_in, multithread = MULTITHREAD, verbose = FALSE, seed = NULL, post_samplename_pattern1 = "_R1.*\\.fastq", post_samplename_pattern2 = "_R2.*\\.fastq", nbases = 1e7){
  if (!is.null(seed)) set.seed(seed)

  fn_fullname <- file.path(dir_in, fn)

  filtFs <- fn_fullname[file.exists(fn_fullname) & grepl(post_samplename_pattern1, fn_fullname)]
  if(length(filtFs) == 0) warning(paste0("runDadaITS: ", "No files found at specified location(s) within ", dir_in, ". Check file path, or post_samplename_pattern argument(s)."))

  # Create outnames
  sample.names <- sapply(strsplit(basename(filtFs), paste0("(",post_samplename_pattern1,")|(",post_samplename_pattern2,")")), `[`, 1)
  names(filtFs) <- sample.names

  # Learn error rates
  errF <- learnErrors(filtFs, nbases=nbases, randomize=TRUE, multithread=multithread, verbose=verbose)

  # Create output vectors
  derepF <- ddF <- vector("list", length(sample.names))
  names(derepF) <- names(ddF) <- sample.names

  # Sample inference and merger of paired-end reads
  for(i in 1:length(sample.names)) {
    sam <- sample.names[[i]]
    cat("Processing:", sam, "\n")
    derepF[[i]] <- derepFastq(filtFs[[i]], verbose=verbose)
    ddF[[i]] <- dada(derepF[[i]], err=errF, multithread=multithread, verbose=verbose)
    cat(paste("Exact sequence variants inferred for sample:", sam,". \n"))
  }
  # rm(derepF);

  # Construct sequence table and remove chimeras
  seqtab <- makeSequenceTable(ddF)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=multithread, verbose=verbose)

  if(verbose){
    cat(paste0("\n\nDimensions of ESV table: ", dim(seqtab)[1], " samples, ", dim(seqtab)[2], " ESVs\n"))
    cat("\nDistribution of sequence lengths:")
    print(table(nchar(getSequences(seqtab))))
    cat(paste0("\n", round(sum(seqtab.nochim)/sum(seqtab), 3)*100, "% reads remain after removal of chimeras"))
  }

  # Track reads through this part of the pipeline
  track <- Reduce(
    function(x, y, ...) transform(merge(x, y, by=0, all = TRUE, ...), row.names=Row.names, Row.names=NULL),
    list(sapply(derepF, getN),
         sapply(ddF, getN),
         rowSums(seqtab.nochim))
  )
  names(track) <- c("derepF", "denoisedF", "nonchim")
  track[is.na(track)] <- 0

  return(list("seqtab" = seqtab,
              "seqtab.nochim" = seqtab.nochim,
              "track" = track))
}


#' Run Dada on R1 ITS Sequences (with metadata)
#'
#' Runs the Dada algorithm to infer sample composition. Currently only supports (R1)
#' sequences. If reverse reads are included, they will be ignored.
#' We use only the R1 reads because attempting to merge with R2
#' reads may result in less accurate representations of the fungal community
#' composition (Pauvert et al., 2019).
#'
#' This implementation is based on Ben Callahan's vignettes at \url{https://benjjneb.github.io/dada2/bigdata.html}
#' and \url{https://benjjneb.github.io/dada2/ITS_workflow.html}.
#'
#' @param fn Full names of input fastq files, including directory. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param meta Metadata downloaded using \code{\link{downloadSequenceMetadata}} that corresponds to the fastq files.
#' @param out_seqtab,out_track (Optional) File locations where sequence table and read-tracking table (respectively) will be saved as csv files. If blank (default), will not save to file.
#' @param remove_chimeras (Optional) Default TRUE. Whether to remove chimeras from the sequence table. Currently only supports removal using \code{\link[dada2]{removeBimeraDenovo}} with the consensus method.
#' @param multithread Default MULTITHREAD in params.R. Whether to use multithreading.
#' @param verbose Default FALSE. Whether to print messages regarding the samples being processed, dimensions of the resulting sequence table, and the distribution of sequence lengths.
#' @param seed (Optional) Integer to use as random seed for reproducibility.
#' @param nbases (Optional) Number of bases to use for learning errors. Default 1e7.
#'
#' @return A list of two elements. \strong{seqtab} is the sequence table, with chimeras removed if remove_chimeras == TRUE, and \strong{track} is a data frame displaying the number of reads remaining for each sample at various points throughout the processing pipeline.
#' @export
#'
#' @examples
#' \dontrun{
#' seqtab.list <- runDada16S(c("sample1_R1.fastq", "sample1_R2.fastq", "sample2_R1.fastq", "sample2_R2.fastq"), meta, seed=1010100)
#' }
runDadaITS2 <- function(fn, meta, out_seqtab = "", out_track = "", remove_chimeras = TRUE, multithread = MULTITHREAD, verbose = FALSE, seed = NULL, nbases=1e7){
  if (!is.null(seed)) set.seed(seed)

  keep_fn <- file.exists(fn) & !duplicated(fn)
  fn <- fn[keep_fn]
  if(length(fn) == 0) warning(paste0("runDada16S: No files found at specified location(s). Check file paths, or input metadata."))

  # Reference metadata to retrieve R1 files
  meta_ext <- matchFastqToMetadata(fn, meta)
  fnFs <- fn[grep("R1", meta_ext$rawDataFileDescription)]

  # Confirm target gene
  if(any(!grepl("16S", meta_ext$targetGene))) warning("You are using runDada16S() on some non-16S files. Did you mean to use runDadaITS()?")

  # Attempt to get DNA sample IDs to use as row names
  dnaSampleIDs <- as.character(meta_ext$dnaSampleID[match(fnFs, meta_ext$file)])
  fn_as_rownames <- any(is.na(dnaSampleIDs))
  if(fn_as_rownames) {
    warning(sum(is.na(dnaSampleIDs)),
            " files do not have associated dnaSampleIDs in the input metadata. ",
            "Therefore, input R1 file names, instead of dnaSampleIDs, will be used as the rownames for the read tracking table.")
  }

  names(fnFs) <- dnaSampleIDs

  # Learn forward and reverse error rates
  errF <- learnErrors(fnFs, nbases=nbases, randomize=TRUE, multithread=multithread, verbose=verbose)

  # Create output vectors
  derepF <- ddF <- vector("list", length(dnaSampleIDs))
  names(derepF) <- names(ddF) <- dnaSampleIDs

  # Sample inference and merger of paired-end reads
  message("runDadaITS progress out of ", length(dnaSampleIDs), " samples:")
  progressbar <- txtProgressBar(min = 0, max = length(dnaSampleIDs), style = 3)
  for(i in 1:length(dnaSampleIDs)) {
    sam <- dnaSampleIDs[[i]]
    if(verbose) cat("Processing:", sam, "\n")

    derepF[[i]] <- derepFastq(fnFs[[i]], verbose=verbose)
    ddF[[i]] <- dada(derepF[[i]], err=errF, multithread=multithread, verbose=verbose)
    if(verbose) cat(paste("Exact sequence variants inferred for sample:", sam,". \n"))
    setTxtProgressBar(progressbar, i)
  }
  close(progressbar)

  # Construct sequence table and remove chimeras
  seqtab <- makeSequenceTable(ddF)
  if(remove_chimeras) seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=multithread, verbose=verbose)


  if(verbose){
    cat(paste0("\nDimensions of ESV table: ", dim(seqtab)[1], " samples, ", dim(seqtab)[2], " ESVs\n"))
    cat("\nDistribution of sequence lengths (should be bimodal for v3v4 region):")
    print(table(nchar(getSequences(seqtab))))
    if(remove_chimeras) cat(paste0("\n", round(sum(seqtab.nochim)/sum(seqtab), 3)*100, "% reads remain after removal of chimeras"))
  }

  # Track reads through this part of the pipeline
  suppressWarnings({
    track <- Reduce(
      function(x, y, ...) transform(merge(x, y, by=0, all = TRUE, ...), row.names=Row.names, Row.names=NULL),
      if(remove_chimeras) {
        list(sapply(derepF, getN),
             sapply(ddF, getN),
             rowSums(seqtab.nochim))
      } else {
        list(sapply(derepF, getN),
             sapply(ddF, getN))
      }
    )
  })
  if(remove_chimeras) {
    names(track) <- c("derepF", "denoisedF", "nonchim")
  } else {
    names(track) <- c("derepF", "denoisedF")
  }
  track[is.na(track)] <- 0

  if(out_seqtab != "") {
    if(!dir.exists(dirname(out_seqtab))) {
      dir.create(dirname(out_seqtab), recursive=TRUE)
      message("Output directory created for sequence table: ", dirname(out_seqtab))
    }
    if(remove_chimeras) {
      write.csv(seqtab.nochim, out_seqtab)
    } else {
      write.csv(seqtab, out_seqtab)
    }
  }
  if(out_track != "") {
    if(!dir.exists(dirname(out_track))) {
      dir.create(dirname(out_track), recursive=TRUE)
      message("Output directory created for read tracking table: ", dirname(out_track))
    }
    write.csv(track, out_track)
  }

  if(remove_chimeras) {
    return(list("seqtab" = seqtab.nochim,
                "track" = track))
  } else {
    return(list("seqtab" = seqtab,
                "track" = track))
  }
}
