# DADA2 workflow for processing NEON ITS raw sequences
# Follows https://benjjneb.github.io/dada2/ITS_workflow.html

# Load parameters from params.R
source("./code/params.R")
path <- file.path(PRESET_OUTDIR_SEQUENCE, "ITS")

# Load libraries
library(dada2)
library(ShortRead)
library(Biostrings)
library(tibble)
library(dplyr)

# Get all run IDs so you can group by them
runs <- regmatches(list.files(path), gregexpr("^run[A-Za-z0-9]*", list.files(path)))
unique_runs <- unique(unlist(runs))

# If SMALL_SUBSET == TRUE, run only the first runID
if (SMALL_SUBSET) {
  loop_length <- 1
} else {
  loop_length <- length(unique_runs)
}

# TODO: Switch from a for-loop to a foreach, or some other parallel process
# system.time({
for (i in c(15)) {
# for (i in 1:loop_length) { # <-- TODO: need to re-run #15 (runBTJKN) and #6 (runBF8M2, which lacks rownames?)
  runID <- unique_runs[i]
  print(paste0("Began processing ", runID, " at ", Sys.time()))
  
  # Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
  fnFs <- sort(list.files(path, pattern=paste0(runID, ".*_R1.fastq"), full.names = TRUE)) # If set to be FALSE, then working directory must contain the files
  fnRs <- sort(list.files(path, pattern=paste0(runID, ".*_R2.fastq"), full.names = TRUE))
  
  # Remove any forward files that don't have reverse counterparts, and vise versa
  # (filterAndTrim will throw an error if fnFs and fnRs have any mismatches)
  basefilenames_Fs <- sub("_R1.fastq","",basename(fnFs))
  basefilenames_Rs <- sub("_R2.fastq","",basename(fnRs))
  rm_from_fnFs <- basefilenames_Fs[which(!(basefilenames_Fs %in% basefilenames_Rs))]
  rm_from_fnRs <- basefilenames_Rs[which(!(basefilenames_Rs %in% basefilenames_Fs))]
  
  for(name in rm_from_fnFs) {
    if(VERBOSE) print(paste(name, "does not have a reverse-reads counterpart. Omitting from this analysis."))
    fnFs <- fnFs[-which(fnFs == paste0(path, "/", name, "_R1.fastq"))]
  }
  for(name in rm_from_fnRs) {
    if(VERBOSE) print(paste(name, "does not have a forward-reads counterpart. Omitting from this analysis."))
    fnRs <- fnRs[-which(fnRs == paste0(path, "/", name, "_R2.fastq"))]
  }
  
  # If SMALL_SUBSET == TRUE,
  # keep only the first two forward-reverse pairs of sequence files
  if(SMALL_SUBSET){
    if(length(fnFs > 2)) fnFs <- fnFs[1:2]
    if(length(fnRs > 2)) fnRs <- fnRs[1:2]
  }
  
  # Get all orientations of primers, for trimming later
  FWD.orients <- allOrients(PRIMER_ITS_FWD)
  REV.orients <- allOrients(PRIMER_ITS_REV)
  
  # “pre-filter” the sequences just to remove those with Ns, but perform no other filtering
  fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
  fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
  # system.time({ # 27.9 s on runB69PP (144 files), 77.5 s on runB69RF (173 files)
    out_filtN <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = MULTITHREAD, compress = FALSE)
  # })
  print(paste0("Finished pre-filtering sequences in ", runID, " at ", Sys.time()))
  
  # This part deviates from the tutorial. Since some samples lose all reads at
  # the pre-filtering stage, it is useful to trim down the samples for the
  # next step: cutadapt
  path.filtN <- file.path(path, "filtN")
  fnFs.filtN2 <- list.files(path.filtN, pattern = paste0(runID, ".*_R1.fastq"), full.names=TRUE)
  fnRs.filtN2 <- list.files(path.filtN, pattern = paste0(runID, ".*_R2.fastq"), full.names=TRUE)
  
  # Count the number of times the forward and reverse primers appear in a
  # set of paired end reads. First file in each sequencing run should
  # be sufficient
  if(VERBOSE) {
    count_primer_orients(fnFs.filtN2[[1]], fnRs.filtN2[[1]], PRIMER_ITS_FWD, PRIMER_ITS_REV)
  }
  # If you see the reverse-complement of the forward primer in the reverse reads (cell [2,4]),
  # or the reverse-complement of the reverse primer in the forward reads (cell [3,4]),
  # it's because the ITS region is short and the read overlaps with the reverse-complement
  # of the other primer.
  
  # TODO: Fix the mixed-orientation reads issue
  
  # Remove primers using cutadapt
  path.cut <- file.path(path, "cutadapt")
  if(!dir.exists(path.cut)) dir.create(path.cut)
  fnFs.cut <- file.path(path.cut, basename(fnFs.filtN2))
  fnRs.cut <- file.path(path.cut, basename(fnRs.filtN2))
  
  FWD.RC <- dada2:::rc(PRIMER_ITS_FWD)
  REV.RC <- dada2:::rc(PRIMER_ITS_REV)
  # Trim FWDPrimer and the reverse-complement of REVPrimer off of R1 (forward reads)
  R1.flags <- paste("-g", PRIMER_ITS_FWD, "-a", REV.RC) 
  # Trim REVPrimer and the reverse-complement of FWDPrimer off of R2 (reverse reads)
  R2.flags <- paste("-G", PRIMER_ITS_REV, "-A", FWD.RC) 
  
  # Run Cutadapt
  # system.time({ # 212.7 s on runB69PP 
  for(i in seq_along(fnFs.filtN2)) { # MAY BE POSSIBLE TO PARALLELIZE THIS STEP
      system2(CUTADAPT_PATH, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                               fnFs.filtN2[i], fnRs.filtN2[i], # input files; fnFs.filtN replaced by fnFs.filtN2, etc.
                               "--minimum-length", "1"), # min length of cutadapted reads: >0 
              stdout = FALSE)
  }
  # })
  print(paste0("Finished stage 1 (of 2) of primer removal using cutadapt in ", runID, " at ", Sys.time()))
  
  # Count primers in first post-cutadapt sample (should all be 0):
  if(VERBOSE) {
    count_primer_orients(fnFs.cut[[1]], fnRs.cut[[1]], PRIMER_ITS_FWD, PRIMER_ITS_REV)
  }
  
  # Since they are not all zero, remove all other orientations of primers
  # EDIT: This is not the best way to handle mixed orientations --
  #       see Issue #5.
  
  path.cut2 <- file.path(path, "cutadapt2")
  if(!dir.exists(path.cut2)) dir.create(path.cut2)
  fnFs.cut2 <- file.path(path.cut2, basename(fnFs.cut))
  fnRs.cut2 <- file.path(path.cut2, basename(fnRs.cut))
  
  # Trim REV and the reverse-complement of FWD off of R1 (forward reads)
  R1.flags.swapped <- paste("-g", PRIMER_ITS_REV, "-a", FWD.RC) 
  # Trim FWD and the reverse-complement of REV off of R2 (reverse reads)
  R2.flags.swapped <- paste("-G", PRIMER_ITS_FWD, "-A", REV.RC) 
  
  # Run Cutadapt
  # system.time({ # 162.6 s on runB69PP
  for(i in seq_along(fnFs.cut)) { # MAY BE POSSIBLE TO PARALLELIZE THIS STEP
    system2(CUTADAPT_PATH, args = c(R1.flags.swapped, R2.flags.swapped, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut2[i], "-p", fnRs.cut2[i], # output files
                               fnFs.cut[i], fnRs.cut[i], # input files; fnFs.filtN replaced by fnFs.filtN2, etc.
                               "--minimum-length", "1"), # min length of cutadapted reads: >0 
            stdout = FALSE)
  }
  # })
  print(paste0("Finished stage 2 (of 2) of primer removal using cutadapt in ", runID, " at ", Sys.time()))
  
  # Check primers
  if(VERBOSE) {
    count_primer_orients(fnFs.cut2[[1]], fnRs.cut2[[1]], PRIMER_ITS_FWD, PRIMER_ITS_REV)
  }
  
  # Forward and reverse fastq filenames have the format:
  cutFs <- sort(list.files(path.cut2, pattern = paste0(runID, ".*_R1.fastq"), full.names = TRUE))
  cutRs <- sort(list.files(path.cut2, pattern = paste0(runID, ".*_R2.fastq"), full.names = TRUE))
  
  # Extract sample names, assuming filenames have format:
  get.sample.name <- function(fname) {
    paste(strsplit(basename(fname), "_")[[1]][1:4], collapse="_") # TODO: this could be made more robust
  }
  sample.names <- unname(sapply(cutFs, get.sample.name))
  if(VERBOSE) head(sample.names)
  
  # Inspect read quality profiles of forward reads #1-2
  if(VERBOSE) plotQualityProfile(cutFs[1:2])
  
  # Inspect read quality profiles of reverse reads #1-2
  if(VERBOSE) plotQualityProfile(cutRs[1:2])
  
  # Filter and trim
  
  # Assigning the filenames for the output of the filtered reads 
  # to be stored as fastq.gz files.
  filtFs <- file.path(path.cut2, "filtered", basename(cutFs))
  filtRs <- file.path(path.cut2, "filtered", basename(cutRs))
  
  # system.time({ # 26.2 s on runB69PP
  out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(MAX_EE_FWD, MAX_EE_REV), 
                       truncQ = TRUNC_Q, minLen = MIN_LEN, compress = TRUE, multithread = MULTITHREAD)  # on windows, set multithread = FALSE
  # })
  print(paste0("Finished filterAndTrim in ", runID, " at ", Sys.time()))

  if(VERBOSE) head(out)
  # TODO: Address the following issue:
  # Why do we lose so many reads? For example, in runB69PP:
  #                                          reads.in reads.out
  # runB69PP_BMI_Plate13WellA12_ITS_R1.fastq     4517       223
  # runB69PP_BMI_Plate13WellA3_ITS_R1.fastq      4513       130
  
  filtFs.out <- list.files(file.path(path.cut2, "filtered"), pattern = paste0(runID, ".*_R1.fastq"), full.names=TRUE)
  filtRs.out <- list.files(file.path(path.cut2, "filtered"), pattern = paste0(runID, ".*_R2.fastq"), full.names=TRUE)
  
  # Learn the error rates
  errF <- learnErrors(filtFs.out, multithread = MULTITHREAD)
  errR <- learnErrors(filtRs.out, multithread = MULTITHREAD)
  print(paste0("Finished learning error rates in ", runID, " at ", Sys.time()))
  
  # Visualize estimated error rates
  if(VERBOSE) plotErrors(errF, nominalQ = TRUE)
  
  # Dereplicate identical reads
  derepFs <- derepFastq(filtFs.out, verbose = TRUE)
  derepRs <- derepFastq(filtRs.out, verbose = TRUE)
  # Name the derep-class objects by the sample names
  sample.names <- unname(sapply(filtFs.out, get.sample.name))
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  print(paste0("Finished dereplication in ", runID, " at ", Sys.time()))
  
  # DADA2's core sample inference algorithm
  dadaFs <- dada(derepFs, err = errF, multithread = MULTITHREAD)
  dadaRs <- dada(derepRs, err = errR, multithread = MULTITHREAD)
  print(paste0("Finished DADA2's core sample inference algorithm in ", runID, " at ", Sys.time()))
  
  # Merge pairs
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
  print(paste0("Finished merging pairs in ", runID, " at ", Sys.time()))
  
  # Construct sequence table
  seqtab <- makeSequenceTable(mergers)
  if(VERBOSE) dim(seqtab)
  
  # Remove chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=MULTITHREAD, verbose=TRUE)
  print(paste0("Finished removing chimeras in ", runID, " at ", Sys.time()))
  
  # Inspect distribution of sequence lengths
  if(VERBOSE) hist(nchar(getSequences(seqtab.nochim)))
  
  # Takes 1299 s (21.7 min) to get to this point from the start of the for loop for runB69PP, with multithreading
  
  # Track reads through pipeline
  getN <- function(x) sum(getUniques(x))
  rownames(out) <- sub("_R1.fastq", "", rownames(out)) # TODO: This could be made more robust
  if(nrow(seqtab.nochim) > 1) {
    track <- cbind(out, 
                   sapply(dadaFs, getN)[rownames(out)], 
                   sapply(dadaRs, getN)[rownames(out)], 
                   sapply(mergers, getN)[rownames(out)], 
                   rowSums(seqtab.nochim)[rownames(out)])
  # If processing a single sample, remove the sapply calls: e.g. replace
  # sapply(dadaFs, getN) with getN(dadaFs)
  } else {
    track <- cbind(out,
                   getN(dadaFs),
                   getN(dadaRs),
                   getN(mergers),
                   rowSums(seqtab.nochim))
  }
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                       "nonchim")
  
  # Save the table which tracks the no. of reads remaining at each stage
  # in the DADA2 pipeline
  path.track <- file.path(path, "track_reads")
  if(!dir.exists(path.track)) dir.create(path.track)
  write.csv(track, file.path(path.track, paste0("track_reads_",runID,".csv")))
  print(paste0("Finished tracking reads through pipeline in ", runID, " at ", Sys.time()))
  
  # Save sequence table associated with this sequencing run
  saveRDS(seqtab.nochim, paste0("./data/NEON_ITS_seqtab_nochim_DL08-13-2019_", runID, ".Rds")) # TODO: May need to include output data file as a parameter in params.R
  print(paste0("Finished saving sequence table of ", runID, " at ", Sys.time()))
  
  # TODO: Need to test the following alternative to full_join:
  # Merge sequence tables
  if(exists("seqtab_joined")) {
    seqtab_joined <- mergeSequenceTables(seqtab_joined, seqtab.nochim)
  } else {
    seqtab_joined <- seqtab.nochim
  }
  
  print(paste0("Finished processing ", runID, " at ", Sys.time()))
}
# })

# Takes 38210 s (10 h 36.8 min) to process runB69RF, runB69RN, runB9994, runBDR3T, and runBF8M2 through the for loop

# Assign taxonomy using the UNITE database
unite.ref <- "./raw_data/tax_ref/sh_general_release_dynamic_02.02.2019.fasta"
taxa_joined <- assignTaxonomy(seqtab_joined, unite.ref, multithread = MULTITHREAD, tryRC = TRUE)
print(paste0("Finished assigning taxonomy in ", runID, " at ", Sys.time()))

taxa.print <- taxa_joined  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
if(VERBOSE) head(taxa.print)

# Saved joined sequence table and joined taxa table
saveRDS(seqtab_joined, paste0("./data/NEON_ITS_seqtab_nochim_DL08-13-2019_", "test", ".Rds")) # TODO: May need to include output data file as a parameter in params.R
saveRDS(taxa_joined, paste0("./data/NEON_ITS_taxa_DL08-13-2019_", "test", ".Rds"))

# Hand off to dada2_to_phyloseq.R
