# DADA2 tutorial
# Follows the ITS workflow https://benjjneb.github.io/dada2/ITS_workflow.html

library(dada2)
library(ShortRead)
library(Biostrings)

path <- "/afs/cats.ucsc.edu/users/b/claraqin/zhulab/NEON_DoB_analysis/data/Illumina/NEON/ITS"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE)) # If set to be FALSE, then working directory must contain the files
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))

# Remove any forward files that don't have reverse counterparts, and vise versa
# (filterAndTrim will throw an error if fnFs and fnRs have any mismatches)
basefilenames_Fs <- sub("_R1.fastq.gz","",basename(fnFs))
basefilenames_Rs <- sub("_R2.fastq.gz","",basename(fnRs))
rm_from_fnFs <- basefilenames_Fs[which(!(basefilenames_Fs %in% basefilenames_Rs))]
rm_from_fnRs <- basefilenames_Rs[which(!(basefilenames_Rs %in% basefilenames_Fs))]
# During this tutorial, I found that "BMI_Plate60WellG10_ITS" exists in the reverse-read
# files but not in the forward-read files. I'll remove this.

for(name in rm_from_fnFs) {
  print(paste(name, "does not have a reverse-reads counterpart. Omitting from this analysis."))
  fnFs <- fnFs[-which(fnFs == paste0(path, "/", name, "_R1.fastq.gz"))]
}
for(name in rm_from_fnRs) {
  print(paste(name, "does not have a forward-reads counterpart. Omitting from this analysis."))
  fnRs <- fnRs[-which(fnRs == paste0(path, "/", name, "_R2.fastq.gz"))]
}


# Identify primers
FWD <- "CTTGGTCATTTAGAGGAAGTAA"  ## CHANGE ME to your forward primer sequence
REV <- "GCTGCGTTCTTCATCGATGC"  ## CHANGE ME...


# Get all orientations of primers, just to be safe
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# “pre-filter” the sequences just to remove those with Ns, but perform no other filtering
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
# "Mismatched forward and reverse sequence files" error occurs many times:
# One of these is at i=1261
# I believe this happens whenever at least one file in the forward-reverse pair
# does not pass the "N" filter.

# Only 4392 out of 11049 ITS fastq files passed?

# This part deviates from the tutorial, but since not all fastq files passed
# the filter, it might make more sense to trim down the filenames for the
# next step: cutadapt
fn2_basenames <- list.files("./data/Illumina/NEON/ITS/filtN", full.names=FALSE)
fnFs2 <- file.path(path, fn2_basenames[grep("_R1.fastq.gz", fn2_basenames)])
fnRs2 <- file.path(path, fn2_basenames[grep("_R2.fastq.gz", fn2_basenames)])
fnFs.filtN2 <- file.path(path, "filtN", fn2_basenames[grep("_R1.fastq.gz", fn2_basenames)])
fnRs.filtN2 <- file.path(path, "filtN", fn2_basenames[grep("_R2.fastq.gz", fn2_basenames)])



# (From tutorial) We are now ready to count the number of times the primers appear in the 
# forward and reverse read, while considering all possible primer orientations. 
# Identifying and counting the primers on one set of paired end FASTQ files is
# sufficient, assuming all the files were created using the same library preparation,
# so we’ll just process the first sample.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN2[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN2[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN2[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN2[[1]]))
# If you see the reverse-complement of the forward primer in the reverse reads (cells [2,4] and [3,4]),
# it's because the ITS region is short and it is reading part of the forward primer.

# Remove primers using cutadapt

cutadapt <- "/afs/cats.ucsc.edu/users/b/claraqin/.local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R


path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs.filtN2))
fnRs.cut <- file.path(path.cut, basename(fnRs.filtN2))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs.filtN2)) {
# for(i in 1:10) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN2[i], fnRs.filtN2[i], # input files; fnFs.filtN replaced by fnFs.filtN2, etc.
                             "--minimum-length", "1")) # min length of cutadapted reads: >0 
}

# Count primers in first post-cutadapt sample (should all be 0):
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Since they are not all zero, remove all other orientations of primers

path.cut2 <- file.path(path, "cutadapt2")
if(!dir.exists(path.cut2)) dir.create(path.cut2)
fnFs.cut2 <- file.path(path.cut2, basename(fnFs.cut))
fnRs.cut2 <- file.path(path.cut2, basename(fnRs.cut))

# Trim REV and the reverse-complement of FWD off of R1 (forward reads)
R1.flags.swapped <- paste("-g", REV, "-a", FWD.RC) 
# Trim FWD and the reverse-complement of REV off of R2 (reverse reads)
R2.flags.swapped <- paste("-G", FWD, "-A", REV.RC) 
# Run Cutadapt
for(i in seq_along(fnFs.cut)) {
# for(i in 1:10) {
  system2(cutadapt, args = c(R1.flags.swapped, R2.flags.swapped, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut2[i], "-p", fnRs.cut2[i], # output files
                             fnFs.cut[i], fnRs.cut[i], # input files; fnFs.filtN replaced by fnFs.filtN2, etc.
                             "--minimum-length", "1")) # min length of cutadapted reads: >0 
}

# Check primers
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut2[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut2[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut2[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut2[[1]]))


# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut2, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut2, pattern = "_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) {
  paste(strsplit(basename(fname), "_")[[1]][1:2], collapse="_")
}
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Inspect read quality profiles of forward reads #1-2
plotQualityProfile(cutFs[1:2])

# Inspect read quality profiles of reverse reads #1-2
plotQualityProfile(cutRs[1:2])

# Filter and trim

# Assigning the filenames for the output of the filtered reads 
# to be stored as fastq.gz files.
filtFs <- file.path(path.cut2, "filtered", basename(cutFs))
filtRs <- file.path(path.cut2, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

filtFs.out <- list.files(paste(path.cut2, "filtered", sep="/"), pattern="_R1.fastq.gz", full.names=TRUE)
filtRs.out <- list.files(paste(path.cut2, "filtered", sep="/"), pattern="_R2.fastq.gz", full.names=TRUE)

# Learn the error rates
errF <- learnErrors(filtFs.out, multithread = TRUE)
errR <- learnErrors(filtRs.out, multithread = TRUE)

# Visualize estimated error rates
plotErrors(errF, nominalQ = TRUE)

# Dereplicate identical reads
derepFs <- derepFastq(filtFs.out, verbose = TRUE)
derepRs <- derepFastq(filtRs.out, verbose = TRUE)
# Name the derep-class objects by the sample names
sample.names <- unname(sapply(filtFs.out, get.sample.name))
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# DADA2's core sample inference algorithm
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# Merge pairs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Inspect distribution of sequence lengths
hist(nchar(getSequences(seqtab.nochim)))

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

### NOTE: unlike the tutorial, we lost a lot of reads at the "filtN" stage,
### which filters out any sequences containing an "N" (ambiguous) base.
### However, the tutorial does acknowledge that a majority of reads may
### be lost in the filtering stage.

# Assign taxonomy using the UNITE database
unite.ref <- "../data/tax_ref/sh_general_release_dynamic_02.02.2019.fasta"
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Hand off to phyloseq



library(phyloseq)
packageVersion("phyloseq")
