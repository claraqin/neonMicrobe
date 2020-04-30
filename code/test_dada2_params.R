# test_dada2_params.R
#
# This script tests the effects of various filterAndTrim and dada2 parameters on
# - proportion / number of reads remaining in a sample
# - diversity estimates of a sample
# - taxonomic resolution of a sample
# 
# Before running this script, you need to have downloaded at least the data
# associated with sequencing run B69PP. You can do this by setting the 'SEQUENCING_RUNS'
# parameter to 'B69PP' in the params.R file, making other adjustments to the parameters
# as needed, and then running 00_new_server_setup.R.

# Load parameters from params.R, tools from utils.R
source("./code/params.R")
source("./code/utils.R")

# Load libraries
library(dada2)
library(ShortRead)
library(Biostrings)
library(tibble)
library(dplyr)
library(vegan)

# Generate filepath names
PATH_ITS <- file.path(PRESET_OUTDIR_SEQUENCE, "ITS")
PATH_UNZIPPED <- file.path(PATH_ITS, "0_unzipped")
PATH_FILTN <- file.path(PATH_ITS, "1_filtN")
PATH_CUT <- file.path(PATH_ITS, "2_cutadapt")
PATH_FILTERED <- file.path(PATH_ITS, "3_filtered")
PATH_SEQTABS <- file.path(PATH_ITS, "4_seqtabs")
PATH_TRACK <- file.path(PATH_ITS, "track_reads")
PATH_TEST <- file.path(PATH_ITS, "test")

# Get all run IDs so you can group by them
unique_runs <- unique(unlist(
  regmatches(list.files(PATH_UNZIPPED), gregexpr("^run[A-Za-z0-9]*", list.files(PATH_UNZIPPED)))
))

##############
# Track proportion of reads remaining after filterAndTrim under different params
i = 18
# This is an arbitrary selection; higher-quality runs include runBTJKN (#15), runC25G9 (#18), runC3CN4 (#19)
runID <- unique_runs[i]
runID
cutFs <- sort(list.files(PATH_CUT, pattern = paste0(runID, ".*_R1.fastq"), full.names = TRUE))
cutRs <- sort(list.files(PATH_CUT, pattern = paste0(runID, ".*_R2.fastq"), full.names = TRUE))

# To cut down on computation time, select 5 samples from the run:
cutFs <- cutFs[round(quantile(1:length(cutFs)))]
cutRs <- cutRs[round(quantile(1:length(cutRs)))]

# Plot quality profiles
for (sample_no in 1:length(cutFs)) {
  gridExtra::grid.arrange(plotQualityProfile(cutFs[[sample_no]]),
                          plotQualityProfile(cutRs[[sample_no]]), nrow=1)  
}

# filterAndTrim using different sets of parameters
param_sets <- c("default","maxEE5","truncQ10","maxEE5_truncQ10")
params <- matrix(c(2,2,2,50,
                   5,5,2,50,
                   2,2,10,50,
                   5,5,10,50), 
                 nrow=4, byrow=TRUE,
                 dimnames = list(param_sets,
                                 c("maxEE.F", "maxEE.R", "truncQ", "minLen")))

out0 <- filterAndTrim(fwd = cutFs, filt = file.path(PATH_TEST, param_sets[1], basename(cutFs)), 
                      rev = cutRs, filt.rev = file.path(PATH_TEST, param_sets[1], basename(cutRs)), 
                      compress = TRUE, multithread = TRUE, maxN = 0,
                      maxEE = params[1,1:2], truncQ = params[1,3], minLen = params[1,4])
out_maxEE5 <- filterAndTrim(cutFs, file.path(PATH_TEST, param_sets[2], basename(cutFs)), cutRs, file.path(PATH_TEST, param_sets[2], basename(cutRs)), compress = TRUE, multithread = TRUE, maxN = 0,
                            maxEE = params[2,1:2], truncQ = params[2,3], minLen = params[2,4])
out_truncQ10 <- filterAndTrim(cutFs, file.path(PATH_TEST, param_sets[3], basename(cutFs)), cutRs, file.path(PATH_TEST, param_sets[3], basename(cutRs)), compress = TRUE, multithread = TRUE, maxN = 0,
                              maxEE = params[3,1:2], truncQ = params[3,3], minLen = params[3,4])
out_maxEE5_truncQ10 <- filterAndTrim(cutFs, file.path(PATH_TEST, param_sets[4], basename(cutFs)), cutRs, file.path(PATH_TEST, param_sets[4], basename(cutRs)), compress = TRUE, multithread = TRUE, maxN = 0,
                                     maxEE = params[4,1:2], truncQ = params[4,3], minLen = params[4,4])

out0_df <- as.data.frame(out0) %>% mutate(prop.out = reads.out / reads.in)
out_maxEE5_df <- as.data.frame(out_maxEE5) %>% mutate(prop.out = reads.out / reads.in)
out_truncQ10_df <- as.data.frame(out_truncQ10) %>% mutate(prop.out = reads.out / reads.in)
out_maxEE5_truncQ10_df <- as.data.frame(out_maxEE5_truncQ10) %>% mutate(prop.out = reads.out / reads.in)

# Number/proportion of reads remaining
out0_df
out_maxEE5_df
out_truncQ10_df
out_maxEE5_truncQ10_df


###############
# Run parts of dada2 (post-filterAndTrim) to see downstream effects on merging, 
# taxonomic resolution, diversity

# Get filenames
filtFs <- list(file.path(PATH_TEST, param_sets[1], basename(cutFs)),
               file.path(PATH_TEST, param_sets[2], basename(cutFs)),
               file.path(PATH_TEST, param_sets[3], basename(cutFs)),
               file.path(PATH_TEST, param_sets[4], basename(cutFs)))
filtFs <- lapply(filtFs, function(x) { x <- x[file.exists(x)] })
filtRs <- list(file.path(PATH_TEST, param_sets[1], basename(cutRs)),
               file.path(PATH_TEST, param_sets[2], basename(cutRs)),
               file.path(PATH_TEST, param_sets[3], basename(cutRs)),
               file.path(PATH_TEST, param_sets[4], basename(cutRs)))
filtRs <- lapply(filtRs, function(x) { x <- x[file.exists(x)] })


n_merged <- list()
prop_merged <- list()
seqtabs <- list()
for(i in 1:length(filtFs)) {
  filtFs.star <- filtFs[[i]]
  filtRs.star <- filtRs[[i]]
  
  set.seed(11001100)
  # Learn the error rates
  errF <- learnErrors(filtFs.star, multithread=MULTITHREAD, nbases = 1e7, randomize=TRUE)
  errR <- learnErrors(filtRs.star, multithread=MULTITHREAD, nbases = 1e7, randomize=TRUE)
  print(paste0("Finished learning error rates in ", param_sets[i], " at ", Sys.time()))
  
  # Visualize estimated error rates
  # if(VERBOSE) plotErrors(errF, nominalQ = TRUE)
  
  # Dereplicate identical reads
  derepFs <- derepFastq(filtFs.star, verbose = TRUE)
  derepRs <- derepFastq(filtRs.star, verbose = TRUE)
  # Name the derep-class objects by the sample names
  sample.names <- unname(sapply(filtFs.star, get.sample.name))
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  print(paste0("Finished dereplication in ", param_sets[i], " at ", Sys.time()))
  
  # DADA2's core sample inference algorithm
  dadaFs <- dada(derepFs, err = errF, multithread = MULTITHREAD)
  dadaRs <- dada(derepRs, err = errR, multithread = MULTITHREAD)
  print(paste0("Finished DADA2's core sample inference algorithm in ", params[i], " at ", Sys.time()))
  
  # Merge pairs
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, returnRejects=TRUE)
  print(paste0("Finished merging pairs in ", param_sets[i], " at ", Sys.time()))
  # rm(derepFs)
  # rm(derepRs)
  
  # If using returnRejects=TRUE in mergePairs(), you can look at the number/proportion
  # of sequences in each sample which successfully merged
  n_merged[[i]] <- unlist(lapply(mergers, function(x) sum(x$accept)))
  prop_merged[[i]] <- unlist(lapply(mergers, function(x) mean(x$accept)))
  
  # Construct sequence table
  # If using returnRejects=TRUE in mergePairs(), you will have to remove the column
  # corresponding to unmerged sequence pairs, "".
  seqtab <- makeSequenceTable(mergers)
  dim(seqtab)
  ind_blank <- which(colnames(seqtab)=="")
  if(length(ind_blank) > 0) {
    seqtab <- seqtab[,-ind_blank]
  }
  
  # Remove chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=MULTITHREAD, verbose=TRUE)
  print(paste0("Finished removing chimeras in ", param_sets[i], " at ", Sys.time()))
  
  # Inspect distribution of sequence lengths
  hist(nchar(getSequences(seqtab.nochim)))
  
  seqtabs[[i]] <- seqtab.nochim
}
names(n_merged) <- param_sets
names(prop_merged) <- param_sets
names(seqtabs) <- param_sets

n_merged
prop_merged


# What's the best way to evaluate effects on diversity estimates?
# Do the parameters affect the shape of the rarefaction curves,
# or just the sequencing depth?
par(mfrow=c(2,2))
lapply(seqtabs, function(x) {
  rarecurve(x, step=50, label=FALSE)
})
par(mfrow=c(1,1))
rarecurve(seqtabs[[4]], sample=50)

# Assign taxonomy using the UNITE database
unite.ref <- "./raw_data/tax_ref/sh_general_release_dynamic_02.02.2019.fasta"
taxas <- list()
for (i in 1:length(seqtabs)) {
  taxas[[i]] <- assignTaxonomy(seqtabs[[i]], unite.ref, multithread = MULTITHREAD, tryRC = TRUE)
}
names(taxas) <- param_sets

taxas.print <- taxas[[1]]
rownames(taxas.print) <- NULL
head(taxas.print)

# What's the best way to evaluate taxonomic resolution?
# Number of non-missing values in the Species column? Not sure
for(i in 1:length(param_sets)) {
  print(paste0("Parameter set: ", param_sets[i]))
  print(paste0("- Nonmissing species-level classifications: ", sum(!is.na(taxas[[i]][,"Species"]))))
  print(paste0("- Total number of ASVs: ", nrow(taxas[[i]])))
}

