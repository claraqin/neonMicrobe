
# Functions for the NEON 16S workflow

trimPrimers16S <- function(fnFs, fnRs, PATH_CUT, PRIMER_16S_FWD, PRIMER_16S_REV, MULTITHREAD#, quiet=T
                           ){
  
  # Create output directory
  if(!dir.exists(PATH_CUT)) dir.create(PATH_CUT)
  
  # Create primer-trimmed file paths in output directory
  fnFs.cut <- file.path(PATH_CUT, basename(fnFs))
  fnRs.cut <- file.path(PATH_CUT, basename(fnRs))
  
  # Just trim the length of fwd and rev primers
  qa.out <- filterAndTrim(fnFs, fnFs.cut, fnRs, fnRs.cut, multithread = MULTITHREAD,
                            matchIDs = TRUE, trimLeft = c(nchar(PRIMER_16S_FWD), nchar(PRIMER_16S_REV)), compress=F)
  
  colnames(qa.out) <- c("input", "trimmed")
  return(qa.out)
}


# Remove any forward files that don't have reverse counterparts, and vise versa
# (filterAndTrim will throw an error if fnFs and fnRs have any mismatches)
remove_unmatched_files <- function(fnFs, fnRs){
  
basefilenames_Fs <- sub("_R1.fastq","",basename(fnFs))
basefilenames_Rs <- sub("_R2.fastq","",basename(fnRs))
rm_from_fnFs <- basefilenames_Fs[which(!(basefilenames_Fs %in% basefilenames_Rs))]
rm_from_fnRs <- basefilenames_Rs[which(!(basefilenames_Rs %in% basefilenames_Fs))]

if (length(c(rm_from_fnFs, rm_from_fnRs)) == 0){ 
  if(VERBOSE == TRUE) cat("All R1/R2 files had counterparts.\n")
  } else {
    for(name in rm_from_fnFs) {
      if(VERBOSE) message(paste(name, "does not have an R2 counterpart. Omitting from this analysis.\n"))
      fnFs <- fnFs[-which(fnFs == paste0(PATH_UNZIPPED, "/", name, "_R1.fastq"))]
    }
    for(name in rm_from_fnRs) {
      if(VERBOSE) message(paste(name, "does not have an R1 counterpart. Omitting from this analysis.\n"))
      fnRs <- fnRs[-which(fnRs == paste0(PATH_UNZIPPED, "/", name, "_R2.fastq"))]
    }
    rm(rm_from_fnFs)
    rm(rm_from_fnRs)
  }
}



qualityFilter16S <- function(PATH_CUT, PATH_FILTERED, MULTITHREAD, MAX_EE_FWD, MAX_EE_REV, TRUNC.LENGTHS = NULL, mean = FALSE){
  fnFs <- sort(list.files(PATH_CUT, pattern = "_R1.fastq.gz|_R1.fastq", full.names = TRUE)) # primer-trimmed 
  fnRs <- sort(list.files(PATH_CUT, pattern = "_R2.fastq.gz|_R2.fastq", full.names = TRUE)) # primer-trimmed 
  
  sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1) # extract sample names
  filtFs <- file.path(PATH_FILTERED, paste0(sample.names, "_filt_R1.fastq.gz")) # create filtered filenames
  filtRs <- file.path(PATH_FILTERED, paste0(sample.names, "_filt_R2.fastq.gz")) # create filtered filenames
  
  n_files <- ifelse(length(unique(sample.names)) > 30, 30, length(unique(sample.names)))
  
  if (is.null(TRUNC.LENGTHS)){
    #### GETTING TRUNCATION LENGTH USING A SUBSET OF QUALITY SCORES ####
    fwd.trunc.lengths <- list()
    rev.trunc.lengths <- list()
  
      fwd.trunc.lengths[1:n_files] <- get_truncation_length(fnFs[1:n_files], verbose = VERBOSE,
                                                                              qscore = 23)
      rev.trunc.lengths[1:n_files] <- get_truncation_length(fnRs[1:n_files], verbose = VERBOSE,
                                                            qscore = 23)
    
    # remove any samples that have low-quality reads early on, to avoid spoiling the whole run
    # idk how to loop this command
    to_remove <- which(unlist(fwd.trunc.lengths) < 100 | unlist(rev.trunc.lengths) < 100)
    if (length(to_remove)>0){
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
    TRUNC.LENGTHS <- c(fwd.trunc.length, rev.trunc.length)
  }
    
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                         multithread=MULTITHREAD,
                         truncLen=TRUNC.LENGTHS,
                         maxEE = c(MAX_EE_FWD, MAX_EE_REV),
                         truncQ=TRUNC_Q, matchIDs = TRUE,
                         compress=TRUE, maxN=0)
    
  return(out)
}




#' Set truncation length
#' 
#' Decides on truncation length for trimmed reads based on quality score means. Default cutoff is a score of 30. 
#' Warns you if the beginning (first 10) bases are low-quality, but returns the first low-quality base after the 10th base.
#' If no bases are below score, returns the last base.
#' 
#'  fl : input files (prints suggested length for each; if used in script, could take the first input, or the lowest)
#'  qscore : default = 30.
#'  n : default = 5e+05 don't know exactly why this matters, but it's part of how 'qa' is called within the dada2 scripts...
#' 
#' 
get_truncation_length <-function (fl, qscore = 30, n = 5e+05, verbose = TRUE){
  trunc_lengths <- data.frame(file = character(0), early_lowqual = numeric(0),
                              trunc_length = numeric(0))
  for (f in fl) {
    srqa <- ShortRead::qa(f, n = n)
    df <- srqa[["perCycle"]]$quality
    means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count, 
                                                          df$Cycle)
    lowqual <- which(means < qscore)
    lowqual_in_first_20 <- length(which(lowqual <= 20))
    lowqual_after_20 <- lowqual[lowqual > 20][1]
    
    trunc_lengths <- rbind(trunc_lengths, data.frame(file = f, 
                                                     early_lowqual = lowqual_in_first_20,
                                                     trunc_length = ifelse(is.na(lowqual_after_20), length(means), lowqual_after_20))) 
    
    if (verbose == TRUE) {
      if(lowqual_in_first_20 > 0){
        cat(paste0("Note: for sample ", basename(f),': ', lowqual_in_first_20, ' base(s) within first 20 bases are below your quality score. Consider trimming left.\n'))
      } 
      if (is.na(lowqual_after_20)){
        cat(paste0(basename(f), ': After first 20 bases, no bases with a mean under your quality score found. Truncate at end of read, base: ', length(means),'\n'))
        
      } else if (!is.na(lowqual_after_20)){
        cat(paste0(basename(f),': After first 20 bases, first mean below your quality score is base: ',lowqual_after_20,'\n'))

      } else "Something's wrong. Inspect function/reads."
    } # end printing to console  
  } # end loop
  return(trunc_lengths$trunc_length)  
} # end function




# Workflow from https://benjjneb.github.io/dada2/bigdata_paired.html
runDada16S <- function(PATH_FILTERED, MULTITHREAD, VERBOSE, seed = NULL, ...){
  if (!is.null(seed)) set.seed(seed)
  
  # File parsing
  filtFs <- list.files(PATH_FILTERED, pattern="_R1.fastq.gz", full.names = TRUE)
  filtRs <- list.files(PATH_FILTERED, pattern="_R2.fastq.gz", full.names = TRUE)
  
  # Keep files with counterpart
  remove_unmatched_files(filtFs, filtRs)

  # Create outnames
  sample.names <- sapply(strsplit(basename(filtFs), "_R"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
  sample.namesR <- sapply(strsplit(basename(filtRs), "_R"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
  if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
  
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  # Learn forward and reverse error rates
  errF <- learnErrors(filtFs, nbases=1e7, randomize=TRUE, multithread=MULTITHREAD)
  errR <- learnErrors(filtRs, nbases=1e7, randomize=TRUE, multithread=MULTITHREAD)
  
  # Create output vectors
  mergers <- vector("list", length(sample.names))
  derepF <- vector("list", length(sample.names))
  derepR <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  names(derepF) <- sample.names
  names(derepR) <- sample.names
  
  # Sample inference and merger of paired-end reads
  for(i in 1:length(sample.names)) {
    sam <- sample.names[[i]]
    cat("Processing:", sam, "\n")
    
    derepF[[i]] <- derepFastq(filtFs[[i]])
    derepR[[i]] <- derepFastq(filtRs[[i]])
    
    ddF <- dada(derepF[[i]], err=errF, multithread=TRUE)
    ddR <- dada(derepR[[i]], err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF[[i]], ddR, derepR[[i]], maxMismatch=1, minOverlap = 6)
    mergers[[i]] <- merger
    cat(paste("Exact sequence variants inferred for sample:",  sam,". \n"))
  }
  #rm(derepF); rm(derepR)
  # Construct sequence table and remove chimeras
  seqtab <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=MULTITHREAD, verbose=VERBOSE)
  
  if(VERBOSE){
  cat(paste0("\n\nDimensions of ESV table: ", dim(seqtab.nochim)[1], " samples, ", dim(seqtab.nochim)[2], " ESVs\n"))
  cat("\nDistribution of sequence lengths (should be bimodal for v3v4 region):")
  print(table(nchar(getSequences(seqtab))))
  cat(paste0("\n", round(sum(seqtab.nochim)/sum(seqtab), 3), " reads remain after removal of chimeras"))
  }
  track <- cbind.data.frame(denoisedF = sapply(derepF, getN),
                            denoisedR = sapply(derepR, getN),
                            merged = sapply(mergers, getN),
                            nonchim = rowSums(seqtab.nochim))
  return(list("seqtab" = seqtab, "seqtab.nochim" = seqtab.nochim, "track" = track))
}

# For tracking reads throughout pipeline
getN <- function(x) sum(getUniques(x))
