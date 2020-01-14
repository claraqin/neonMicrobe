# DADA2 workflow for processing NEON ITS raw sequences
# Follows https://benjjneb.github.io/dada2/tutorial.html

library(dada2)
library(ShortRead)
library(Biostrings)

path <- "/data/ZHULAB/NEON_DOB/Illumina/NEON/16S"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE)) # If set to be FALSE, then working directory must contain the files
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

# Remove any forward files that don't have reverse counterparts, and vise versa
# (filterAndTrim will throw an error if fnFs and fnRs have any mismatches)
basefilenames_Fs <- sub("_R1.fastq","",basename(fnFs))
basefilenames_Rs <- sub("_R2.fastq","",basename(fnRs))
rm_from_fnFs <- basefilenames_Fs[which(!(basefilenames_Fs %in% basefilenames_Rs))]
rm_from_fnRs <- basefilenames_Rs[which(!(basefilenames_Rs %in% basefilenames_Fs))]

for(name in rm_from_fnFs) {
  print(paste(name, "does not have a reverse-reads counterpart. Omitting from this analysis."))
  fnFs <- fnFs[-which(fnFs == paste0(path, "/", name, "_R1.fastq"))]
}
for(name in rm_from_fnRs) {
  print(paste(name, "does not have a forward-reads counterpart. Omitting from this analysis."))
  fnRs <- fnRs[-which(fnRs == paste0(path, "/", name, "_R2.fastq"))]
}

sample.names <- sub("_R1.fastq","",basename(fnFs))

# Inspect read quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Filter and trim

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

# Error in filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(240, 160),  : 
# These are the errors (up to 5) encountered in individual cores...
# Error in (function (fn, fout, maxN = c(0, 0), truncQ = c(2, 2), truncLen = c(0,  : 
# Mismatched forward and reverse sequence files: 49729, 28061.
# Error in (function (fn, fout, maxN = c(0, 0), truncQ = c(2, 2), truncLen = c(0,  : 
# Mismatched forward and reverse sequence files: 49729, 28061.
# Error in (function (fn, fout, maxN = c(0, 0), truncQ = c(2, 2), truncLen = c(0,  : 
# Mismatched forward and reverse sequence files: 49729, 28061.
# Error in (function (fn, fout, maxN = c(0, 0), truncQ = c(2, 2), truncLen = c(0,  : 
# Mismatched forward and reverse sequence files: 49729, 28061.
# Error in (function (fn, fout, maxN = c(0, 0), truncQ = c(2, 2), truncLen = c(0,  : 
# Mismatched forward and reverse sequence files: 49729, 28061.
# In addition: Warning message:
# In mclapply(seq_len(n), do_one, mc.preschedule = mc.preschedule,  :
# scheduled cores 22 encountered errors in user code, all values of the jobs will be affected

# "Mismatched forward and reverse sequence files" error occurs many times:
# I believe this happens whenever at least one file in the forward-reverse pair
# does not pass the "N" filter.

# 11250 out of 11378 16S fastq files passed

# This part deviates from the tutorial, but since not all fastq files passed
# the filter, it might make more sense to trim down the filenames for the
# next step
filt2_basenames <- list.files(file.path(path,"filtered"), full.names=FALSE)
filtFs2 <- file.path(path, "filtered", filt2_basenames[grep("_F_filt.fastq", filt2_basenames)])
filtRs2 <- file.path(path, "filtered", filt2_basenames[grep("_R_filt.fastq", filt2_basenames)])


# Learn the error rates
errF <- learnErrors(filtFs2, multithread=TRUE)
errR <- learnErrors(filtRs2, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)


# DADA2's core sample inference algorithm
dadaFs <- dada(filtFs2, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs2, err=errR, multithread=TRUE)
dadaFs[[1]]


# Merge pairs
mergers <- mergePairs(dadaFs, filtFs2, dadaRs, filtRs2, verbose=TRUE)
head(mergers[[1]])


# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# What proportion remains after removal of chimeras?
sum(seqtab.nochim)/sum(seqtab)


# # Track reads through pipeline
# getN <- function(x) sum(getUniques(x))
# track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
# colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
# rownames(track) <- sample.names
# head(track)


# Assign taxonomy using the Silva reference database
taxa <- assignTaxonomy(seqtab.nochim, "./raw_data/tax_ref/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Save OTU table and taxonomic table as RDS files
# to hand off to dada2_to_phyloseq.R
saveRDS(seqtab.nochim, "./data/NEON_16S_seqtab_nochim_DL08-13-2019.Rds")
saveRDS(taxa, "./data/NEON_16S_taxa_DL08-13-2019.Rds")
