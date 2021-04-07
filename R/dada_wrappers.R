
#' Validate Input/Output Arguments to DADA2 Wrappers
#'
#' Internal function for validating the input and output path arguments
#' to any of the DADA2 wrapper functions. Keeps only existing input
#' files. Output directory is modifiable by processing batches via
#' \code{\link{NEONMICROBE_DIR_OUTPUTS}}.
#'
#' @param fn Base names of input fastq files. If inputs are not base names (i.e. if they include directory paths), the directory paths will be removed. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param in_subdir Subdirectory name from which to retrieve input fastq files. Enter "raw" for raw sequence files, or any other character string to specify a subdirectory within \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S. To specify a directory outside \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S, use the 'in_explicitdir' argument.
#' @param out_subdir Subdirectory name where output fastq files will be written. Enter any character string to specify a subdirectory within \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S. If the directory does not exist, it will be created. To specify a directory outside \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S, use the 'out_explicitdir' argument.
#' @param in_explicitdir,out_explicitdir Directory names to use instead of 'in_subdir' and 'out_subdir', if static directory names or directories outside of \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S are desired. Not recommended for use within processing batches.
#' @param targetGene Either "16S" or "ITS".
#'
#' @return List.
#' - fn: Full paths to input files.
#' - fn_out: Full paths where each file will be written after processing.
#' - in_dir: Full path to input directory.
#' - out_dir: Full path to output directory.
validateIO <- function(fn, in_subdir, out_subdir, in_explicitdir, out_explicitdir, targetGene = c("16S", "ITS"),
                       validate_input_only = FALSE) {
  targetGene <- match.arg(targetGene)

  # Validate input filenames
  if(!all(dirname(fn) %in% c(".", ""))) {
    warning("At least one input fastq filename was not a base name (i.e. contained a directory path). ",
            "Directory paths in input filenames will be ignored.")
  }

  # Validate input/output dirs
  if((!missing(in_subdir) && !is.null(in_explicitdir)) ||
     (!missing(out_subdir) && !is.null(out_explicitdir))) {
    warning("Both implicit and explicit directory names were passed to function. Explicit directory ",
            "names will take precedence.")
  }
  if(is.null(in_explicitdir)) {
    if(identical(in_subdir, "raw")) { # special directory path if in_subdir=="raw"
      in_dir <- file.path(NEONMICROBE_DIR_SEQUENCE(), targetGene)
    } else {
      in_dir <- file.path(NEONMICROBE_DIR_MIDPROCESS(), targetGene, in_subdir)
    }
  } else {
    in_dir <- in_explicitdir
  }
  if(!identical(validate_input_only, TRUE)) {
    if(is.null(out_explicitdir)) {
      out_dir <- file.path(NEONMICROBE_DIR_MIDPROCESS(), targetGene, out_subdir)
    } else {
      out_dir <- out_explicitdir
    }
  }

  # Get input files
  fn <- file.path(in_dir, basename(fn))
  keep_fn <- file.exists(fn) & !duplicated(fn)
  fn <- fn[keep_fn]
  if(length(fn) == 0) {
    warning(paste0("No fastq files found at specified location(s). Check file paths."))
    return(invisible(NULL))
  }

  # Create output directory and output filenames
  if(!identical(validate_input_only, TRUE)) {
    dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
    fn_out <- file.path(out_dir, basename(fn))
  }

  if(!identical(validate_input_only, TRUE)) {
    return(list(
      fn = fn,
      fn_out = fn_out,
      in_dir = in_dir,
      out_dir = out_dir
    ))
  } else {
    return(list(
      fn = fn,
      in_dir = in_dir
    ))
  }
}

#' Trim Primers from 16S Sequences (with metadata)
#'
#' Trims primers from 16S sequences using \code{\link[dada2]{filterAndTrim}}.
#' This function assumes that each read begins with its full primer sequence
#' and operates by truncating the beginning of each read by the length of its
#' primer sequence.
#'
#' @param fn Base names of input fastq files. If inputs are not base names (i.e. if they include directory paths), the directory paths will be removed. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param in_subdir Subdirectory name from which to retrieve input fastq files. Enter "raw" for raw sequence files, or any other character string to specify a subdirectory within \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S. To specify a directory outside \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S, use the 'in_explicitdir' argument.
#' @param out_subdir Subdirectory name where output fastq files will be written. Enter any character string to specify a subdirectory within \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S. If the directory does not exist, it will be created. To specify a directory outside \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S, use the 'out_explicitdir' argument.
#' @param meta The output of \code{\link{downloadSequenceMetadata}}. Must be provided as either the data.frame returned by \code{\link{downloadSequenceMetadata}} or as a filepath to the csv file produced by \code{\link{downloadSequenceMetadata}}.
#' @param in_explicitdir,out_explicitdir Directory names to use instead of 'in_subdir' and 'out_subdir', if static directory names or directories outside of \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S are desired. Not recommended for use within processing batches.
#' @param primer_16S_fwd,primer_16S_rev DNA sequences of 16S forward and reverse primer, respectively
#' @param multithread Default FALSE. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#'
#' @return (Invisibly) Integer matrix denoting the number of reads remaining after primer-trimming for each input file.
#' @export
#'
#' @examples
trimPrimers16S <- function(fn, in_subdir, out_subdir, meta, in_explicitdir = NULL, out_explicitdir = NULL, primer_16S_fwd = "CCTACGGGNBGCASCAG", primer_16S_rev = "GACTACNVGGGTATCTAATCC", multithread = FALSE) {

  # Validate sequence metadata
  checkArgsAgainstBatchParams(meta = "SEQMETA_FILE")
  meta <- readSequenceMetadata(meta)

  # Validate input and output arguments
  validated_args <- validateIO(fn, in_subdir, out_subdir, in_explicitdir, out_explicitdir, targetGene="16S")
  fn <- validated_args$fn
  fn_out <- validated_args$fn_out
  in_dir <- validated_args$in_dir
  out_dir <- validated_args$out_dir
  if(length(fn) == 0) {
    return(invisible(NULL))
  }

  # Get metadata matching files
  meta_ext <- matchFastqToMetadata(fn, meta)

  # Reference metadata to retrieve R1 and R2 files
  fn_pairs <- getPairedFastqFiles(fn, meta_ext[,-1], value=FALSE)
  fnFs <- fn[fn_pairs[[1]]]
  fnRs <- fn[fn_pairs[[2]]]

  # Confirm target gene
  if(any(!grepl("16S", meta_ext$targetGene))) warning("You are using trimPrimers16S() on some non-16S files. Did you mean to use trimPrimersITS()?")

  # Attempt to get DNA sample IDs to use as row names
  dnaSampleIDs <- as.character(meta_ext$dnaSampleID[fn_pairs[[1]]])
  fn_as_rownames <- any(is.na(dnaSampleIDs))
  if(fn_as_rownames) {
    warning(sum(is.na(dnaSampleIDs)),
            " files do not have associated dnaSampleIDs in the input metadata. ",
            "Therefore, input R1 file names, instead of dnaSampleIDs, will be used as the rownames for the read tracking table.")
  }

  # Confirm output filenames
  trimFs <- fn_out[fn_pairs[[1]]]
  trimRs <- fn_out[fn_pairs[[2]]]

  # Just trim the length of fwd and rev primers
  qa.out <- dada2::filterAndTrim(fnFs, trimFs, fnRs, trimRs, multithread = multithread,
                                 matchIDs = TRUE, trimLeft = c(nchar(primer_16S_fwd), nchar(primer_16S_rev)),
                                 compress=TRUE)

  message("Input fastq files trimmed and written to ", out_dir)

  colnames(qa.out) <- c("input", "trimmed")
  if(fn_as_rownames) {
    rownames(qa.out) <- fnFs
  } else {
    rownames(qa.out) <- dnaSampleIDs
  }
  return(invisible(qa.out))
}


#' Trim Primers from ITS Sequences (with metadata)
#'
#' Trims primers from ITS sequences using cutadapt. Cutadapt must be installed in order for this to work. Currently only supports forward-read (R1) sequences. If reverse reads are included, they will be ignored.
#'
#' @param fn Base names of input fastq files. If inputs are not base names (i.e. if they include directory paths), the directory paths will be removed. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param in_subdir Subdirectory name from which to retrieve input fastq files. Enter "raw" for raw sequence files, or any other character string to specify a subdirectory within \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/ITS. To specify a directory outside \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/ITS, use the 'in_explicitdir' argument.
#' @param out_subdir Subdirectory name where output fastq files will be written. Enter any character string to specify a subdirectory within \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/ITS. If the directory does not exist, it will be created. To specify a directory outside \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/ITS, use the 'out_explicitdir' argument.
#' @param meta The output of \code{\link{downloadSequenceMetadata}}. Must be provided as either the data.frame returned by \code{\link{downloadSequenceMetadata}} or as a filepath to the csv file produced by \code{\link{downloadSequenceMetadata}}.
#' @param cutadapt_path Default CUTADAPT_PATH in params.R. Path to cutadapt on your file system.
#' @param in_explicitdir,out_explicitdir Directory names to use instead of 'in_subdir' and 'out_subdir', if static directory names or directories outside of \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/ITS are desired. Not recommended for use within processing batches.
#' @param primer_ITS_fwd,primer_ITS_rev DNA sequence of the ITS forward and reverse primer, respectively.
#' @param very_verbose Default FALSE. Whether to print output from cutadapt.
#' @param discard_untrimmed Default FALSE. Whether to discard reads where a primer could not be found, leaving only those reads in which a primer has been trimmed.
#'
#' @return No value is returned.
#' @export
#'
#' @examples
trimPrimersITS <- function(fn, in_subdir, out_subdir, meta, cutadapt_path, in_explicitdir = NULL, out_explicitdir = NULL, primer_ITS_fwd = "CTTGGTCATTTAGAGGAAGTAA", primer_ITS_rev = "GCTGCGTTCTTCATCGATGC", very_verbose=FALSE, discard_untrimmed=FALSE) {

  # Validate sequence metadata
  checkArgsAgainstBatchParams(meta = "SEQMETA_FILE")
  meta <- readSequenceMetadata(meta)

  # Validate input and output arguments
  validated_args <- validateIO(fn, in_subdir, out_subdir, in_explicitdir, out_explicitdir, targetGene="ITS")
  fn <- validated_args$fn
  fn_out <- validated_args$fn_out
  in_dir <- validated_args$in_dir
  out_dir <- validated_args$out_dir
  if(length(fn) == 0) {
    return(invisible(NULL))
  }

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
  trimFs <- fn_out[grep("R1", meta_ext$rawDataFileDescription)]

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

  message("Input fastq files trimmed and written to ", out_dir)
}




#' Filter 16S Sequences (with metadata)
#'
#' Applies a quality filter to 16S sequence fastq files via the \code{\link[dada2]{filterAndTrim}} function.
#' It is assumed that both forward- and reverse-read files are included.
#'
#' @param fn Base names of input fastq files. If inputs are not base names (i.e. if they include directory paths), the directory paths will be removed. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param in_subdir Subdirectory name from which to retrieve input fastq files. Enter "raw" for raw sequence files, or any other character string to specify a subdirectory within \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S. To specify a directory outside \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S, use the 'in_explicitdir' argument.
#' @param out_subdir Subdirectory name where output fastq files will be written. Enter any character string to specify a subdirectory within \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S. If the directory does not exist, it will be created. To specify a directory outside \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S, use the 'out_explicitdir' argument.
#' @param meta The output of \code{\link{downloadSequenceMetadata}}. Must be provided as either the data.frame returned by \code{\link{downloadSequenceMetadata}} or as a filepath to the csv file produced by \code{\link{downloadSequenceMetadata}}.
#' @param in_explicitdir,out_explicitdir Directory names to use instead of 'in_subdir' and 'out_subdir', if static directory names or directories outside of \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S are desired. Not recommended for use within processing batches.
#' @param multithread Default FALSE. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#' @param verbose Default FALSE. Whether to print messages associated with automated selection of truncation length (truncLen).
#' @param trunc_qscore Default 23. If truncLen is not provided, specifies the mean quality score at which point to truncate each read. The behavior of this argument has not been extensively tested, so it is generally recommended that truncLen be used instead.
#' @param ... Other arguments to be passed to \code{\link[dada2]{filterAndTrim}}, such as maxEE, truncLen, and minLen. By default, uses \code{\link[dada2]{filterAndTrim}}'s default values. See documentation for \code{\link[dada2]{filterAndTrim}} for more details.
#'
#' @seealso \code{\link[dada2]{filterAndTrim}}
#'
#' @return (Invisibly) Two-column matrix displaying the number of reads in input vs. output for each file.
#' @export
qualityFilter16S <- function(fn, in_subdir, out_subdir, meta, in_explicitdir = NULL, out_explicitdir = NULL,
                             multithread = FALSE, verbose = FALSE, trunc_qscore = 23, ...){

  # Validate sequence metadata
  checkArgsAgainstBatchParams(meta = "SEQMETA_FILE")
  meta <- readSequenceMetadata(meta)

  # Validate input and output arguments
  validated_args <- validateIO(fn, in_subdir, out_subdir, in_explicitdir, out_explicitdir, targetGene="16S")
  fn <- validated_args$fn
  fn_out <- validated_args$fn_out
  in_dir <- validated_args$in_dir
  out_dir <- validated_args$out_dir
  if(length(fn) == 0) {
    return(invisible(NULL))
  }

  # Get metadata matching files
  meta_ext <- matchFastqToMetadata(fn, meta)

  # Reference metadata to retrieve R1 and R2 files
  fn_pairs <- getPairedFastqFiles(fn, meta_ext[,-1], value=FALSE)
  fnFs <- fn[fn_pairs[[1]]]
  fnRs <- fn[fn_pairs[[2]]]

  # Confirm target gene
  if(any(!grepl("16S", meta_ext$targetGene))) warning("You are using qualityFilter16S() on some non-16S files. Did you mean to use qualityFilterITS()?")

  # Attempt to get DNA sample IDs to use as row names
  dnaSampleIDs <- as.character(meta_ext$dnaSampleID[fn_pairs[[1]]])
  fn_as_rownames <- any(is.na(dnaSampleIDs))
  if(fn_as_rownames) {
    warning(sum(is.na(dnaSampleIDs)),
            " files do not have associated dnaSampleIDs in the input metadata. ",
            "Therefore, input R1 file names, instead of dnaSampleIDs, will be used as the rownames for the read tracking table.")
  }

  # Confirm output filenames
  filtFs <- fn_out[fn_pairs[[1]]]
  filtRs <- fn_out[fn_pairs[[2]]]

  # Read in dots args
  dots <- list(...)

  # If truncLen is NOT in user-provided arguments
  if (!("truncLen" %in% names(dots))) {
    # GETTING TRUNCATION LENGTH USING A SUBSET OF QUALITY SCORES
    n_files <- min(length(fnFs), 30)

    fwd.trunc.lengths <- list()
    rev.trunc.lengths <- list()
    fwd.trunc.lengths[1:n_files] <- getTruncationLength(fnFs[1:n_files], verbose = verbose,
                                                        qscore = trunc_qscore)
    rev.trunc.lengths[1:n_files] <- getTruncationLength(fnRs[1:n_files], verbose = verbose,
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
  out <- do.call(dada2::filterAndTrim, arguments)

  message("Input fastq files filtered and written to ", out_dir)

  if(fn_as_rownames) {
    rownames(out) <- fnFs
  } else {
    rownames(out) <- dnaSampleIDs
  }

  return(invisible(out))
}




#' Filter ITS Sequences (with metadata)
#'
#' Applies a quality filter to ITS sequence fastq files via the \code{\link[dada2]{filterAndTrim}} function.
#' Currently only supports filtering forward-read (R1) sequences. If reverse reads are included, they will be ignored.
#'
#' @param fn Base names of input fastq files. If inputs are not base names (i.e. if they include directory paths), the directory paths will be removed. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param in_subdir Subdirectory name from which to retrieve input fastq files. Enter "raw" for raw sequence files, or any other character string to specify a subdirectory within \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S. To specify a directory outside \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S, use the 'in_explicitdir' argument.
#' @param out_subdir Subdirectory name where output fastq files will be written. Enter any character string to specify a subdirectory within \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S. If the directory does not exist, it will be created. To specify a directory outside \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S, use the 'out_explicitdir' argument.
#' @param meta The output of \code{\link{downloadSequenceMetadata}}. Must be provided as either the data.frame returned by \code{\link{downloadSequenceMetadata}} or as a filepath to the csv file produced by \code{\link{downloadSequenceMetadata}}.
#' @param in_explicitdir,out_explicitdir Directory names to use instead of 'in_subdir' and 'out_subdir', if static directory names or directories outside of \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S are desired. Not recommended for use within processing batches.
#' @param multithread Default FALSE. Whether to use multithreading. Note that Windows does not support multithreading in this function because it uses mclapply, so this argument must be set to FALSE on Windows systems.
#' @param ... Other arguments to be passed to \code{\link[dada2]{filterAndTrim}}, such as maxEE, truncLen, and minLen. By default, uses \code{\link[dada2]{filterAndTrim}}'s default values. See documentation for \code{\link[dada2]{filterAndTrim}} for more details.
#'
#' @return (Invisibly) Two-column matrix displaying the number of reads in input vs. output for each file.
#' @export
qualityFilterITS <- function(fn, in_subdir, out_subdir, meta, in_explicitdir = NULL, out_explicitdir = NULL,
                             multithread = FALSE, ...){
  # Validate sequence metadata
  checkArgsAgainstBatchParams(meta = "SEQMETA_FILE")
  meta <- readSequenceMetadata(meta)

  # Validate input and output arguments
  validated_args <- validateIO(fn, in_subdir, out_subdir, in_explicitdir, out_explicitdir, targetGene="ITS")
  fn <- validated_args$fn
  fn_out <- validated_args$fn_out
  in_dir <- validated_args$in_dir
  out_dir <- validated_args$out_dir
  if(length(fn) == 0) {
    return(invisible(NULL))
  }

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
  filtFs <- fn_out[grep("R1", meta_ext$rawDataFileDescription)]

  # Run quality filter
  out <- dada2::filterAndTrim(fnFs, filtFs,
                              compress = TRUE,
                              multithread = multithread,
                              ...)

  message("Input fastq files filtered and written to ", out_dir)

  if(fn_as_rownames) {
    rownames(out) <- fnFs
  } else {
    rownames(out) <- dnaSampleIDs
  }
  return(invisible(out))
}




#' Run Dada on Paired-End 16S Sequences (with metadata)
#'
#' Runs the Dada algorithm to infer sample composition from paired-end 16S fastq files.
#' This implementation is based on Ben Callahan's vignette at \url{https://benjjneb.github.io/dada2/bigdata_paired.html}.
#'
#' @param fn Base names of input fastq files. If inputs are not base names (i.e. if they include directory paths), the directory paths will be removed. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param in_subdir Subdirectory name from which to retrieve input fastq files. Enter "raw" for raw sequence files, or any other character string to specify a subdirectory within \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S. To specify a directory outside \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S, use the 'in_explicitdir' argument.
#' @param meta The output of \code{\link{downloadSequenceMetadata}}. Must be provided as either the data.frame returned by \code{\link{downloadSequenceMetadata}} or as a filepath to the csv file produced by \code{\link{downloadSequenceMetadata}}.
#' @param out_seqtab,out_track File locations where copies of the sequence table and read-tracking table (respectively) will be written. By default (NULL), these are \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S/3_seqtabs/asv_16s_[timestamp].Rds and \code{\link{NEONMICROBE_DIR_TRACKREADS}}/16S/dada_16s_[timestamp].csv. If no copy should be saved, set the corresponding argument to FALSE.
#' @param in_explicitdir Directory name to use instead of 'in_subdir', if static directory name or directory outside of \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S is desired. Not recommended for use within processing batches.
#' @param remove_chimeras (Optional) Default TRUE. Whether to remove chimeras from the sequence table. Currently only supports removal using \code{\link[dada2]{removeBimeraDenovo}} with the consensus method.
#' @param multithread Default FALSE. Whether to use multithreading.
#' @param verbose Default FALSE. Whether to print messages regarding the samples being processed, dimensions of the resulting sequence table, and the distribution of sequence lengths.
#' @param seed (Optional) Integer to use as random seed for reproducibility.
#' @param nbases (Optional) Number of bases to use for learning errors. Default 1e7.
#'
#' @return (Invisibly) A list of two elements. \strong{seqtab} is the sequence table, with chimeras removed if remove_chimeras == TRUE, and \strong{track} is a data frame displaying the number of reads remaining for each sample at various points throughout the processing pipeline.
#' @export
#'
#' @examples
runDada16S <- function(fn, in_subdir, meta, out_seqtab = NULL, out_track = NULL, in_explicitdir = NULL,
                       remove_chimeras = TRUE, multithread = FALSE, verbose = FALSE, seed = NULL, nbases = 1e7) {

  # Validate sequence metadata
  checkArgsAgainstBatchParams(meta = "SEQMETA_FILE")
  meta <- readSequenceMetadata(meta)

  # Validate input arguments
  validated_args <- validateIO(fn, in_subdir=in_subdir, in_explicitdir=in_explicitdir, targetGene="16S",
                               validate_input_only=TRUE)
  fn <- validated_args$fn
  if(length(fn) == 0) {
    return(invisible(NULL))
  }

  # Validate output arguments
  if(is.null(out_seqtab)) {
    out_seqtab <- file.path(NEONMICROBE_DIR_MIDPROCESS(), "16S", "3_seqtabs",
                            paste0("asv_16s_", sub(" ", "_", gsub(":", "", Sys.time())), ".Rds"))
  }
  if(is.null(out_track)) {
    out_track <- file.path(NEONMICROBE_DIR_TRACKREADS(), "16S",
                            paste0("dada_16s_", sub(" ", "_", gsub(":", "", Sys.time())), ".csv"))
  }

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

  # Set seed
  if(!is.null(seed)) set.seed(seed)

  # Learn forward and reverse error rates
  errF <- dada2::learnErrors(fnFs, nbases=nbases, randomize=TRUE, multithread=multithread, verbose=verbose)
  errR <- dada2::learnErrors(fnRs, nbases=nbases, randomize=TRUE, multithread=multithread, verbose=verbose)

  # Create output vectors
  derepF <- derepR <- ddF <- ddR <- mergers <- vector("list", length(dnaSampleIDs))
  names(derepF) <- names(derepR) <- names(ddF) <- names(ddR) <- names(mergers) <- dnaSampleIDs

  # Sample inference and merger of paired-end reads
  message("runDada16S progress out of ", length(dnaSampleIDs), " samples:")
  progressbar <- txtProgressBar(min = 0, max = length(dnaSampleIDs), style = 3)
  for(i in 1:length(dnaSampleIDs)) {
    sam <- dnaSampleIDs[[i]]
    if(verbose) cat("Processing:", sam, "\n")

    derepF[[i]] <- dada2::derepFastq(fnFs[[i]])
    derepR[[i]] <- dada2::derepFastq(fnRs[[i]])

    ddF[[i]] <- dada2::dada(derepF[[i]], err=errF, multithread=multithread, verbose=verbose)
    ddR[[i]] <- dada2::dada(derepR[[i]], err=errR, multithread=multithread, verbose=verbose)
    if(verbose) {
      mergers[[i]] <- dada2::mergePairs(ddF[[i]], derepF[[i]], ddR[[i]], derepR[[i]],
                                        maxMismatch=1, minOverlap = 6, verbose=TRUE)
    } else {
      mergers[[i]] <- suppressMessages(dada2::mergePairs(ddF[[i]], derepF[[i]], ddR[[i]],
                                                         derepR[[i]], maxMismatch=1,
                                                         minOverlap = 6, verbose=FALSE))
    }
    if(verbose) cat(paste("Exact sequence variants inferred for sample:", sam,". \n"))
    setTxtProgressBar(progressbar, i)
  }
  close(progressbar)

  # Construct sequence table and remove chimeras
  if(verbose) {
    seqtab <- dada2::makeSequenceTable(mergers)
  } else {
    seqtab <- suppressMessages(dada2::makeSequenceTable(mergers))
  }
  if(remove_chimeras) {
    seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus",
                                               multithread=multithread,
                                               verbose=verbose)
  }

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

  if(!identical(out_seqtab, FALSE)) {
    if(!dir.exists(dirname(out_seqtab))) {
      dir.create(dirname(out_seqtab), recursive=TRUE)
      message("Output directory created for sequence table: ", dirname(out_seqtab))
    }
    if(remove_chimeras) {
      saveRDS(seqtab.nochim, out_seqtab)
    } else {
      saveRDS(seqtab, out_seqtab)
    }
    message("Sequence table saved to ", out_seqtab)
  }
  if(!identical(out_track, FALSE)) {
    if(!dir.exists(dirname(out_track))) {
      dir.create(dirname(out_track), recursive=TRUE)
      message("Output directory created for read tracking table: ", dirname(out_track))
    }
    write.csv(track, out_track)
    message("Read tracking table saved to ", out_track)
  }

  if(remove_chimeras) {
    return(invisible(
      list("seqtab" = seqtab.nochim,
           "track" = track)))
  } else {
    return(invisible(
      list("seqtab" = seqtab,
           "track" = track)))
  }
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
#' @param fn Base names of input fastq files. If inputs are not base names (i.e. if they include directory paths), the directory paths will be removed. Files that do not exist will be ignored; however, if all files do not exist, this function will issue a warning.
#' @param in_subdir Subdirectory name from which to retrieve input fastq files. Enter "raw" for raw sequence files, or any other character string to specify a subdirectory within \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S. To specify a directory outside \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S, use the 'in_explicitdir' argument.
#' @param meta The output of \code{\link{downloadSequenceMetadata}}. Must be provided as either the data.frame returned by \code{\link{downloadSequenceMetadata}} or as a filepath to the csv file produced by \code{\link{downloadSequenceMetadata}}.
#' @param out_seqtab,out_track File locations where copies of the sequence table and read-tracking table (respectively) will be written. By default (NULL), these are \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S/3_seqtabs/asv_16s_[timestamp].Rds and \code{\link{NEONMICROBE_DIR_TRACKREADS}}/16S/dada_16s_[timestamp].csv. If no copy should be saved, set the corresponding argument to FALSE.
#' @param in_explicitdir Directory name to use instead of 'in_subdir', if static directory name or directory outside of \code{\link{NEONMICROBE_DIR_MIDPROCESS}}/16S is desired. Not recommended for use within processing batches.
#' @param remove_chimeras (Optional) Default TRUE. Whether to remove chimeras from the sequence table. Currently only supports removal using \code{\link[dada2]{removeBimeraDenovo}} with the consensus method.
#' @param multithread Default FALSE. Whether to use multithreading.
#' @param verbose Default FALSE. Whether to print messages regarding the samples being processed, dimensions of the resulting sequence table, and the distribution of sequence lengths.
#' @param seed (Optional) Integer to use as random seed for reproducibility.
#' @param nbases (Optional) Number of bases to use for learning errors. Default 1e7.
#'
#' @return (Invisibly) A list of two elements. \strong{seqtab} is the sequence table, with chimeras removed if remove_chimeras == TRUE, and \strong{track} is a data frame displaying the number of reads remaining for each sample at various points throughout the processing pipeline.
#' @export
#'
#' @examples
runDadaITS <- function(fn, in_subdir, meta, out_seqtab = NULL, out_track = NULL, in_explicitdir = NULL,
                       remove_chimeras = TRUE, multithread = FALSE, verbose = FALSE, seed = NULL, nbases = 1e7) {

  # Validate sequence metadata
  checkArgsAgainstBatchParams(meta = "SEQMETA_FILE")
  meta <- readSequenceMetadata(meta)

  # Validate input arguments
  validated_args <- validateIO(fn, in_subdir=in_subdir, in_explicitdir=in_explicitdir, targetGene="ITS",
                               validate_input_only=TRUE)
  fn <- validated_args$fn
  if(length(fn) == 0) {
    return(invisible(NULL))
  }

  # Validate output arguments
  if(is.null(out_seqtab)) {
    out_seqtab <- file.path(NEONMICROBE_DIR_MIDPROCESS(), "16S", "4_seqtabs",
                            paste0("asv_its_", sub(" ", "_", gsub(":", "", Sys.time())), ".Rds"))
  }
  if(is.null(out_track)) {
    out_track <- file.path(NEONMICROBE_DIR_TRACKREADS(), "16S",
                           paste0("dada_its_", sub(" ", "_", gsub(":", "", Sys.time())), ".csv"))
  }

  # Reference metadata to retrieve R1 files
  meta_ext <- matchFastqToMetadata(fn, meta)
  fnFs <- fn[grep("R1", meta_ext$rawDataFileDescription)]

  # Confirm target gene
  if(any(!grepl("16S", meta_ext$targetGene))) warning("You are using runDadaITS() on some non-ITS files. Did you mean to use runDada16S()?")

  # Attempt to get DNA sample IDs to use as row names
  dnaSampleIDs <- as.character(meta_ext$dnaSampleID[match(fnFs, meta_ext$file)])
  fn_as_rownames <- any(is.na(dnaSampleIDs))
  if(fn_as_rownames) {
    warning(sum(is.na(dnaSampleIDs)),
            " files do not have associated dnaSampleIDs in the input metadata. ",
            "Therefore, input R1 file names, instead of dnaSampleIDs, will be used as the rownames for the read tracking table.")
  }

  names(fnFs) <- dnaSampleIDs

  # Set seed
  if(!is.null(seed)) set.seed(seed)

  # Learn forward and reverse error rates
  errF <- dada2::learnErrors(fnFs, nbases=nbases, randomize=TRUE, multithread=multithread, verbose=verbose)

  # Create output vectors
  derepF <- ddF <- vector("list", length(dnaSampleIDs))
  names(derepF) <- names(ddF) <- dnaSampleIDs

  # Sample inference and merger of paired-end reads
  message("runDadaITS progress out of ", length(dnaSampleIDs), " samples:")
  progressbar <- txtProgressBar(min = 0, max = length(dnaSampleIDs), style = 3)
  for(i in 1:length(dnaSampleIDs)) {
    sam <- dnaSampleIDs[[i]]
    if(verbose) cat("Processing:", sam, "\n")

    derepF[[i]] <- dada2::derepFastq(fnFs[[i]], verbose=verbose)
    ddF[[i]] <- dada2::dada(derepF[[i]], err=errF, multithread=multithread, verbose=verbose)
    if(verbose) cat(paste("Exact sequence variants inferred for sample:", sam,". \n"))
    setTxtProgressBar(progressbar, i)
  }
  close(progressbar)

  # Construct sequence table and remove chimeras
  seqtab <- dada2::makeSequenceTable(ddF)
  if(remove_chimeras) {
    seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus",
                                               multithread=multithread,
                                               verbose=verbose)
  }


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

  if(!identical(out_seqtab, FALSE)) {
    if(!dir.exists(dirname(out_seqtab))) {
      dir.create(dirname(out_seqtab), recursive=TRUE)
      message("Output directory created for sequence table: ", dirname(out_seqtab))
    }
    if(remove_chimeras) {
      saveRDS(seqtab.nochim, out_seqtab)
    } else {
      saveRDS(seqtab, out_seqtab)
    }
    message("Sequence table saved to ", out_seqtab)
  }
  if(!identical(out_track, FALSE)) {
    if(!dir.exists(dirname(out_track))) {
      dir.create(dirname(out_track), recursive=TRUE)
      message("Output directory created for read tracking table: ", dirname(out_track))
    }
    write.csv(track, out_track)
    message("Read tracking table saved to ", out_track)
  }

  if(remove_chimeras) {
    return(invisible(
      list("seqtab" = seqtab.nochim,
           "track" = track)))
  } else {
    return(invisible(
      list("seqtab" = seqtab,
           "track" = track)))
  }
}


