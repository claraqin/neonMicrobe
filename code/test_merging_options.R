# test_merging_options.R
#
# This script tests the impact of different paired-end read-merging
# options on the number of remaining reads in the output sequence table.
# Inspired by https://github.com/benjjneb/dada2/issues/279

# Merging options:
# 1. R1 only
# 2. Merge R1 and R2 where possible; drop unmerged reads
# 3. Just concatenate R1 and R2
# 4. (TODO:) Merge R1 and R2 where possible; concatenate where merge fails

# NOTE: One of the issues with the current setup is that for the R1-only case,
# I also run filterAndTrim on the R1 reads only. This means that there are more
# reads going into the DADA denoising algorithm for the R1-only case than the
# other cases. I should separate this and start with a common set of samples
# for the DADA algorithm.

# Before running this script, you need to have run 00_new_server_setup.R

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
library(phyloseq)
library(ggplot2)

# Generate filepath names
PATH_ITS <- file.path(PRESET_OUTDIR_SEQUENCE, "ITS")
PATH_UNZIPPED <- file.path(PATH_ITS, "0_unzipped")
PATH_FILTN <- file.path(PATH_ITS, "1_filtN")
PATH_CUT <- file.path(PATH_ITS, "2_cutadapt")
PATH_FILTERED <- file.path(PATH_ITS, "3_filtered")
PATH_SEQTABS <- file.path(PATH_ITS, "4_seqtabs")
PATH_TRACK <- file.path(PATH_ITS, "track_reads")
PATH_TEST <- file.path(PATH_ITS, "test")

runIDs <- c("B69PP", "C25G9")
merging_opts <- c("r1_only", "merge_and_drop", "just_concat")
# merging_opts <- c("r1_only", "merge_and_drop", "just_concat", "merge_and_concat")

##############
# Prep work for the merging option "r1_only"

# Runs to use: (arbitrary selection, though C25G9 is higher quality, B69PP is lower quality)
cutFs <- sort(list.files(PATH_CUT, pattern = paste0("(",paste0("(", runIDs, ")", collapse="|"),")", ".*_R1.fastq"), full.names = TRUE))
cutRs <- sort(list.files(PATH_CUT, pattern = paste0("(",paste0("(", runIDs, ")", collapse="|"),")", ".*_R2.fastq"), full.names = TRUE))

# To cut down on computation time, select up to 100 samples from the runs, up to 100/length(runIDs) from each run.
if(length(cutFs) > 100) {
  cutFs_subset <- c()
  cutRs_subset <- c()
  for(i in 1:length(runIDs)) {
    cutFs_runID <- cutFs[grep(runIDs[i], cutFs)]
    cutRs_runID <- cutRs[grep(runIDs[i], cutRs)]
    cutFs_subset <- c(cutFs_subset, cutFs_runID[round(quantile(1:length(cutFs_runID), probs=seq(0,1,length.out=100/length(runIDs))))])
    cutRs_subset <- c(cutRs_subset, cutRs_runID[round(quantile(1:length(cutRs_runID), probs=seq(0,1,length.out=100/length(runIDs))))])
  }
  cutFs <- cutFs_subset
  cutRs <- cutRs_subset
}

# Plot quality profiles
gridExtra::grid.arrange(plotQualityProfile(cutFs, aggregate=TRUE),
                        plotQualityProfile(cutRs, aggregate=TRUE), nrow=1)

# Create output filenames for the filterAndTrim outputs where only R1
# reads will be used vs. where R1 and R2 will be filtered in conjunction
# with each other
filtFs_r1 <- file.path(PATH_TEST, "merging_opt_r1_only", basename(cutFs))
filtFs_paired <- file.path(PATH_TEST, "merging_opt_paired", basename(cutFs))
filtRs_paired <- file.path(PATH_TEST, "merging_opt_paired", basename(cutRs))

# Run filterAndTrim on samples with either r1 only or paired reads
# FOR THIS ANALYSIS, WE SIMPLY USE THE DEFAULT PARAMETER VALUES FROM params.R
out_list <- list()
out_list[[1]] <- filterAndTrim(
  fwd = cutFs, filt = filtFs_r1,
  compress = TRUE, multithread = TRUE, maxN = 0,
  maxEE = MAX_EE_FWD, truncQ = TRUNC_Q, minLen = MIN_LEN
)
out_list[[2]] <- filterAndTrim(
  fwd = cutFs, filt = filtFs_paired, rev = cutRs, filt.rev = filtRs_paired,
  compress = TRUE, multithread = TRUE, maxN = 0,
  maxEE = c(MAX_EE_FWD, MAX_EE_REV), truncQ = TRUNC_Q, minLen = MIN_LEN
)
out_list[[1]] <- out_list[[1]] %>% as.data.frame() %>% mutate(prop.out = reads.out / reads.in)
out_list[[2]] <- out_list[[2]] %>% as.data.frame() %>% mutate(prop.out = reads.out / reads.in)
names(out_list) <- c("r1_only", "paired")

saveRDS(out_list, "./data/merging_opts_filterAndTrim_out_list.Rds")
out_list <- readRDS("./data/merging_opts_filterAndTrim_out_list.Rds")


#####################
# Run DADA2

# This section will consist of two levels of nested for-loops
# - i = 1:2, iterates over sequencing runs
# - j = either 1:2, iterates over whether R1-only or paired-read cases,
#       or 1:4, iterates over merging options

# And the output objects will be nested lists with corresponding structures:
dadaFs_list <- list(list(), list()) # One sub-list for each of two sequencing runs, each element corresponds to the R1-only vs. paired-read cases
dadaRs_list <- list(list(), list()) # One sub-list for each of two sequencing runs, each element corresponds to the R1-only vs. paired-read cases
mergers_list <- list(list(), list()) # One sub-list for each of two sequencing runs, each element corresponds to a merging option
n_merged <- list(list(), list()) # One sub-list for each of two sequencing runs, each element corresponds to a merging option
prop_merged <- list(list(), list()) # One sub-list for each of two sequencing runs, each element corresponds to a merging option
seqtabs <- list(list(), list()) # One sub-list for each of two sequencing runs, each element corresponds to a merging option

t1 <- Sys.time()
for(i in 1:length(runIDs)) {

  # (1/2)
  # First handle the R1-only case
  filtFs_r1_selrunID <- filtFs_r1[grep(runIDs[i], filtFs_r1)] # select files matching runID
  filtFs_r1_selrunID <- filtFs_r1_selrunID[file.exists(filtFs_r1_selrunID)] # make sure the files exist
  set.seed(11001100)
  errF <- learnErrors(filtFs_r1_selrunID, multithread=MULTITHREAD, nbases = 1e7, randomize=TRUE) # learn error rates
  derepFs <- derepFastq(filtFs_r1_selrunID, verbose = TRUE) # dereplicate
  sample.names <- unname(sapply(filtFs_r1_selrunID, get.sample.name))
  names(derepFs) <- sample.names
  dadaFs <- dada(derepFs, err = errF, multithread = MULTITHREAD) # DADA algorithm
  dadaFs_list[[i]][[1]] <- dadaFs
  seqtab <- makeSequenceTable(derepFs) # construct sequence table
  rm(derepFs)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=MULTITHREAD, verbose=TRUE) # remove chimeras
  seqtabs[[i]][[1]] <- seqtab.nochim

  # (2/2)
  # Then handle the paired-read cases
  filtFs_paired_selrunID <- filtFs_paired[grep(runIDs[i], filtFs_paired)] # select files matching runID
  filtRs_paired_selrunID <- filtRs_paired[grep(runIDs[i], filtRs_paired)]
  filtFs_paired_selrunID <- filtFs_paired_selrunID[file.exists(filtFs_paired_selrunID)] # make sure the files exist
  filtRs_paired_selrunID <- filtRs_paired_selrunID[file.exists(filtRs_paired_selrunID)]
  set.seed(11001100)
  errF <- learnErrors(filtFs_paired_selrunID, multithread=MULTITHREAD, nbases = 1e7, randomize=TRUE) # learn error rates
  errR <- learnErrors(filtRs_paired_selrunID, multithread=MULTITHREAD, nbases = 1e7, randomize=TRUE)
  derepFs <- derepFastq(filtFs_paired_selrunID, verbose = TRUE) # dereplicate
  derepRs <- derepFastq(filtRs_paired_selrunID, verbose = TRUE)
  sample.names <- unname(sapply(filtFs_paired_selrunID, get.sample.name))
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  dadaFs <- dada(derepFs, err = errF, multithread = MULTITHREAD) # DADA algorithm
  dadaRs <- dada(derepRs, err = errR, multithread = MULTITHREAD)
  dadaFs_list[[i]][[2]] <- dadaFs
  dadaRs_list[[i]][[2]] <- dadaRs
  mergers_list[[i]][[2]] <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, returnRejects=TRUE)
  mergers_list[[i]][[3]] <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, returnRejects=TRUE, justConcatenate=TRUE)
  rm(derepFs)
  rm(derepRs)

  for(j in 2:length(merging_opts)) {
    # If using returnRejects=TRUE in mergePairs(), you can look at the number/proportion
    # of sequences in each sample which successfully merged
    n_merged[[i]][[j]] <- unlist(lapply(mergers_list[[i]][[j]], function(x) sum(x$accept)))
    prop_merged[[i]][[j]] <- unlist(lapply(mergers_list[[i]][[j]], function(x) mean(x$accept)))

    # Construct sequence table
    # If using returnRejects=TRUE in mergePairs(), you will have to remove the column
    # corresponding to unmerged sequence pairs, "".
    seqtab <- makeSequenceTable(mergers_list[[i]][[j]])
    ind_blank <- which(colnames(seqtab)=="")
    if(length(ind_blank) > 0) {
      seqtab <- seqtab[,-ind_blank]
    }

    # Remove chimeras
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=MULTITHREAD, verbose=TRUE)

    # Inspect distribution of sequence lengths
    # hist(nchar(getSequences(seqtab.nochim)))

    seqtabs[[i]][[j]] <- seqtab.nochim
  }
}
t2 <- Sys.time() # took 18.5 hours
names(dadaFs_list[[1]]) <- c("r1_only", "paired")
names(dadaFs_list[[2]]) <- c("r1_only", "paired")
names(dadaRs_list[[1]]) <- c("r1_only", "paired")
names(dadaRs_list[[2]]) <- c("r1_only", "paired")
names(mergers_list[[1]]) <- merging_opts
names(mergers_list[[2]]) <- merging_opts
names(n_merged[[1]]) <- merging_opts
names(n_merged[[2]]) <- merging_opts
names(prop_merged[[1]]) <- merging_opts
names(prop_merged[[2]]) <- merging_opts
names(seqtabs[[1]]) <- merging_opts
names(seqtabs[[2]]) <- merging_opts


saveRDS(dadaFs_list, "./data/merging_opts_dadaFs_list.Rds")
saveRDS(dadaRs_list, "./data/merging_opts_dadaRs_list.Rds")
saveRDS(mergers_list, "./data/merging_opts_mergers_list.Rds")
saveRDS(seqtabs, "./data/merging_opts_seqtabs_list.Rds")
saveRDS(n_merged, "./data/merging_opts_n_merged.Rds")
saveRDS(prop_merged, "./data/merging_opts_prop_merged.Rds")



# n_merged # This is the number of ASVs (as partitioned by dada) that successfully merged
# prop_merged # This is proportion of ASVs (as partitioned by dada) that successfully merged
# # The below is the proportion of forward reads that were assigned an ASV
# prop_Fs_mapped_to_asv <- list(list(), list())
# prop_Fs_mapped_to_asv[[1]] <- lapply(dadaFs_list[[1]], function(x) lapply(x, function(y) mean(!is.na(y$map))))
# prop_Fs_mapped_to_asv[[2]] <- lapply(dadaFs_list[[2]], function(x) lapply(x, function(y) mean(!is.na(y$map))))
# prop_Fs_mapped_to_asv_mat <- matrix(
#   unlist(prop_Fs_mapped_to_asv),
#   ncol=100,nrow=16,dimnames=list(param_sets, basename(cutFs))
# )
# head(prop_Fs_mapped_to_asv_mat)
#
# # Rate of successfully merged ASVs per pre-filterAndTrim read
# # (Not sure if this is meaningful -- treats denoising algorithm as a black box)
# merged_variants_per_input_read <- lapply(n_merged, function(x) x/out_list[[1]][,"reads.in"])
# merged_variants_per_input_read



# Assign taxonomy using the UNITE database
# Warning: this can take several days. If available, load from file

# First rename each seqtab to reflect the merging option that created it
seqtabs_renamed <- seqtabs
for(i in 1:length(seqtabs_renamed)) {
  for(j in 1:length(seqtabs_renamed[[i]])) {
    rownames(seqtabs_renamed[[i]][[j]]) <- paste0(merging_opts[j], "_", rownames(seqtabs[[i]][[j]]))
  }
}
seqtab_joinedRuns <- list() # Meaning that the runs have been joined together.
                            # First element corresponds to merging and dropping unsuccessful merges,
                            # second element corresponds to just concatenating
for(i in 1:length(runIDs)) {
  seqtab_joinedRuns[[i]] <- mergeSequenceTables(tables=list(seqtabs_renamed[[i]][[1]], seqtabs_renamed[[i]][[2]], seqtabs_renamed[[i]][[3]]))
}

# Now assign taxonomy
t3 <- Sys.time()
taxas_joinedRuns <- list()
for (i in 1:length(seqtab_joinedRuns)) {
  taxas_joinedRuns[[i]] <- assignTaxonomy(seqtab_joinedRuns[[i]], UNITE_DB_PATH, multithread = MULTITHREAD, tryRC = TRUE)
}
t4 <- Sys.time() # TAXONOMY ASSIGNMENT TOOK 5.888 DAYS (for two sequence tables)

# taxas.print <- taxas_joinedRuns[[1]]
# rownames(taxas.print) <- NULL
# head(taxas.print)

saveRDS(taxas_joinedRuns, "./data/merging_opts_taxas_joinedRuns.Rds")

taxas_joinedRuns <- readRDS("./data/merging_opts_taxas_joinedRuns.Rds")
