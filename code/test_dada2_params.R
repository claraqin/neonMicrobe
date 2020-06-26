# test_dada2_params.R
#
# This script tests the effects of various filterAndTrim and dada2 parameters on
# - proportion / number of reads remaining in a sample
# - diversity estimates of a sample
# - taxonomic resolution of a sample

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

##############
# Runs to use: (arbitrary selection, though C25G9 is higher quality, B69PP is lower quality)
runIDs <- c("B69PP", "C25G9")
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

# Generalize previous version's setup to allow testing filterAndTrim on
# multiple sets of parameters at a time. The following tests 16 sets at a time:
grid <- expand.grid(2^seq(1,4), 2^seq(1,4))
params <- matrix(
  c(grid[,1],     # maxEE.F
    grid[,1],     # maxEE.R
    grid[,2]),     # truncQ
  byrow=FALSE, ncol=3,
  dimnames = list(NULL, c("maxEE.F", "maxEE.R", "truncQ"))
)
params

param_sets <- apply(params, 1, function(x) paste(c(rbind(c("maxEE.F", "maxEE.R", "truncQ"), x)), collapse="_"))

# Create output filenames:
# Since there are 16 parameter sets, and 100 files in the sample,
# each list has 16 list elements, and each list element consists of 100 filenames.
filtFs <- lapply(param_sets, function(x) file.path(PATH_TEST, x, basename(cutFs)))
filtRs <- lapply(param_sets, function(x) file.path(PATH_TEST, x, basename(cutRs)))

# Run filterAndTrim on samples with all 16 parameter sets
out_list <- list()
system.time({
for (i in 1:nrow(params)) {
  out <- filterAndTrim(
    fwd = cutFs, filt = filtFs[[i]], rev = cutRs, filt.rev = filtRs[[i]], 
    compress = TRUE, multithread = TRUE, maxN = 0,
    maxEE = params[i,1:2], truncQ = params[i,3], minLen = 50
  )
  out_list[[i]] <- as.data.frame(out) %>% mutate(prop.out = reads.out / reads.in)
}
})

names(out_list) <- param_sets

# Code for plotting moved to test_dada2_params_plots.Rmd

# These plots don't tell us much, because there is no guarantee that 
# any of these reads are of sufficient quality to be assigned taxonomy,
# or that they will merge.

###############
# Run parts of dada2 (post-filterAndTrim) to see downstream effects on merging, 
# taxonomic resolution, diversity

# Ensure that all of the output filenames exist
filtFs <- lapply(filtFs, function(x) { x <- x[file.exists(x)] })
filtRs <- lapply(filtRs, function(x) { x <- x[file.exists(x)] })

n_merged1 <- list() # For run #1
n_merged2 <- list() # For run #2
n_merged <- list(n_merged1, n_merged2)
prop_merged1 <- list() # For run #1
prop_merged2 <- list() # For run #2
prop_merged <- list(prop_merged1, prop_merged2)
seqtabs1 <- list() # For run #1
seqtabs2 <- list() # For run #2
seqtabs <- list(seqtabs1, seqtabs2)
dadaFs_list1 <- list() # For run #1, forward
dadaFs_list2 <- list() # For run #2, forward
dadaFs_list <- list(dadaFs_list1, dadaFs_list2)
dadaRs_list1 <- list() # For run #1, reverse
dadaRs_list2 <- list() # For run #2, reverse
dadaRs_list <- list(dadaRs_list1, dadaRs_list2)
system.time({
for(i in 1:length(filtFs)) { 
  for(j in 1:length(runIDs)) { 
  # Retrieve only those files associated with the appropriate parameter set and runID
  filtFs.star <- filtFs[[i]][grep(runIDs[j], filtFs[[i]])]
  filtRs.star <- filtRs[[i]][grep(runIDs[j], filtRs[[i]])]
  
  set.seed(11001100)
  # Learn the error rates
  errF <- learnErrors(filtFs.star, multithread=MULTITHREAD, nbases = 1e7, randomize=TRUE)
  errR <- learnErrors(filtRs.star, multithread=MULTITHREAD, nbases = 1e7, randomize=TRUE)
  print(paste0("Finished learning error rates in ", param_sets[i], runID[j], " at ", Sys.time()))
  
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
  
  ## TODO: Record intermediate metric here: the number of
  ##       sequence variants partitioned from the full set
  ##       of reads after the denoising algorithm.
  dadaFs_list[[j]][[i]] <- dadaFs
  dadaRs_list[[j]][[i]] <- dadaRs
  
  # Merge pairs
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, returnRejects=TRUE)
  print(paste0("Finished merging pairs in ", param_sets[i], " at ", Sys.time()))
  # rm(derepFs)
  # rm(derepRs)

  # If using returnRejects=TRUE in mergePairs(), you can look at the number/proportion
  # of sequences in each sample which successfully merged
  n_merged[[j]][[i]] <- unlist(lapply(mergers, function(x) sum(x$accept)))
  prop_merged[[j]][[i]] <- unlist(lapply(mergers, function(x) mean(x$accept)))

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
  print(paste0("Finished removing chimeras in ", param_sets[i], runID[j], "at", Sys.time()))

  # Inspect distribution of sequence lengths
  # hist(nchar(getSequences(seqtab.nochim)))

  seqtabs[[j]][[i]] <- seqtab.nochim
}
}
})
names(n_merged[[1]]) <- param_sets
names(n_merged[[2]]) <- param_sets
names(prop_merged[[1]]) <- param_sets
names(prop_merged[[2]]) <- param_sets
names(dadaFs_list[[1]]) <- param_sets
names(dadaFs_list[[2]]) <- param_sets
names(dadaRs_list[[1]]) <- param_sets
names(dadaRs_list[[2]]) <- param_sets
names(seqtabs[[1]]) <- param_sets
names(seqtabs[[2]]) <- param_sets

n_merged # This is the number of ASVs (as partitioned by dada) that successfully merged
prop_merged # This is proportion of ASVs (as partitioned by dada) that successfully merged
# The below is the proportion of forward reads that were assigned an ASV
prop_Fs_mapped_to_asv <- list(list(), list()) 
prop_Fs_mapped_to_asv[[1]] <- lapply(dadaFs_list[[1]], function(x) lapply(x, function(y) mean(!is.na(y$map))))
prop_Fs_mapped_to_asv[[2]] <- lapply(dadaFs_list[[2]], function(x) lapply(x, function(y) mean(!is.na(y$map))))
prop_Fs_mapped_to_asv_mat <- matrix(
  unlist(prop_Fs_mapped_to_asv),
  ncol=100,nrow=16,dimnames=list(param_sets, basename(cutFs))
)
head(prop_Fs_mapped_to_asv_mat)

# Rate of successfully merged ASVs per pre-filterAndTrim read
# (Not sure if this is meaningful -- treats denoising algorithm as a black box)
merged_variants_per_input_read <- lapply(n_merged, function(x) x/out_list[[1]][,"reads.in"])
merged_variants_per_input_read

# Assign taxonomy using the UNITE database
# Warning: this can take several days. If available, load from file
unite.ref <- "./raw_data/tax_ref/sh_general_release_dynamic_02.02.2019.fasta"
taxas <- list()
for (i in 1:length(seqtabs)) {
  taxas[[i]] <- assignTaxonomy(seqtabs[[i]], unite.ref, multithread = MULTITHREAD, tryRC = TRUE)
}
names(taxas) <- param_sets

taxas.print <- taxas[[1]]
rownames(taxas.print) <- NULL
head(taxas.print)

# Save the data thus far
saveRDS(cutFs, "./data/sensitivity_cutFs.Rds")
saveRDS(cutRs, "./data/sensitivity_cutRs.Rds")
saveRDS(params, "./data/sensitivity_params.Rds")
saveRDS(out_list, "./data/sensitivity_filterAndTrim_out_list.Rds")
saveRDS(prop_Fs_mapped_to_asv, "./data/sensitivity_prop_Fs_mapped_seprun.Rds")
saveRDS(dadaFs_list, "./data/sensitivity_dadaFs_list_seprun.Rds")
saveRDS(dadaRs_list, "./data/sensitivity_dadaRs_list_seprun.Rds")
saveRDS(seqtabs, "./data/sensitivity_seqtabs_list_seprun.Rds")
saveRDS(taxas, "./data/sensitivity_taxas.Rds")
saveRDS(n_merged, "./data/sensitivity_n_merged_seprun.Rds")
saveRDS(prop_merged, "./data/sensitivity_prop_merged_seprun.Rds")

# (Moved plotting code to test_dada2_params_plots.Rmd)

# These results indicate that the best parameters for ensuring
# a high number of successfully merged ASVs are low values
# for truncQ (2-8) and maxEE (2-8).


###################
# Evaluate effects on taxonomic resolution.
# What's the best way to evaluate taxonomic resolution?
# Number of taxa without a species-level assignment?
physeqs <- list() # Create a list of physeq objects
abund_spp_classification <- list() # Total abundance of ASVs with species-level classifications
prop_spp_classification <- list() # Proportion-by-abundance of ASVs with species-level classification
n_spp_classification <- list() # Number of unique ASVs with species-level classification
for(i in 1:length(seqtabs)) {
  physeqs[[i]] <- phyloseq(otu_table(seqtabs[[i]],taxa_are_rows=FALSE), tax_table(taxas[[i]]))
  abund_spp_classification[[i]] <- sum(otu_table(subset_taxa(physeqs[[i]], !is.na(Species))))
  prop_spp_classification[[i]] <- abund_spp_classification[[i]] / sum(otu_table(physeqs[[i]]))
  n_spp_classification[[i]] <- sum(!is.na(tax_table(physeqs[[i]])[,"Species"]))
}

# Remake summary_df with additional columns
summary_df <- 
  data.frame(
    maxEE = params[,1],
    truncQ = params[,3],
    # prop.Fs.mapped = unlist(lapply(prop_Fs_mapped_to_asv, function(x) median(unlist(x)))),
    # prop.merged = unlist(lapply(prop_merged, function(x) median(x))),
    # rate.merged = unlist(lapply(merged_variants_per_input_read, function(x) median(x))),
    abund.spp.classification = unlist(abund_spp_classification),
    prop.spp.classification = unlist(prop_spp_classification),
    n.spp.classification = unlist(n_spp_classification)
  )

# (Moved plotting code to test_dada2_params_plots.Rmd)

# FINDING: This is intuitive; as filtering parameters become more
# stringent, the proportion of ASVs which are assigned to the species
# level increases. However, this does not mean that the absolute
# number of spp-identified ASVs increases; only the spp identification 
# rate.

##################
# What's the best way to evaluate effects on diversity estimates?
# Do the parameters affect the shape of the rarefaction curves,
# or just the sequencing depth?
par(mfrow=c(2,2))
lapply(seqtabs, function(x) {
  rarecurve(x, step=50, label=FALSE)
})
par(mfrow=c(1,1))
rarecurve(seqtabs[[4]], sample=50)

# Use phyloseq's estimate_richness function?
?estimate_richness
estimate_richness(physeqs[[1]])
# Ben Callahan recommends not using diversity estimates that rely on singletons.
# This leaves Shannon and Simpson's indices. https://github.com/benjjneb/dada2/issues/214
estimate_richness(physeqs[[1]], measures=c("Observed","Shannon"))
estimate_richness(physeqs[[1]], measures=c("Observed","Shannon"), split=FALSE) # aggregate into one measure

obsrich_list <- list()
shannon_list <- list()
for(i in 1:length(physeqs)) {
  div <- suppressWarnings(
    estimate_richness(physeqs[[i]], measures=c("Observed","Shannon"), split=FALSE)
  )
  obsrich_list[[i]] <- div[1,1]
  shannon_list[[i]] <- div[1,2]
}

# (Moved code for plotting to test_dada2_params_plots.Rmd)

##############
# Evaluate effects on beta-diversity inference, 
# using ordinations and permANOVA

# (Some of this code is redundant with test_dada2_params_plots.Rmd)

# First rename each seqtab to reflect the params that created it
seqtabs_renamed <- seqtabs
for(i in 1:length(seqtabs_renamed)) {
  rownames(seqtabs_renamed[[i]]) <- paste0(param_sets[i], "_", rownames(seqtabs[[i]]))
}
# Then combine all seqtabs into ONE physeq
seqtab_joined <- mergeSequenceTables(tables=seqtabs_renamed)
rm(seqtabs_renamed)

# Sample data (parameters)
paramset_index <- unlist(lapply(seq_along(seqtabs), function(i) { rep(i, nrow(seqtabs[[i]]))}))
sampledata <- as.data.frame(params[paramset_index,]) %>%
  mutate(sample = unname(unlist(lapply(seqtabs, rownames))))
rownames(sampledata) <- rownames(seqtab_joined)

# Taxa table
# Warning: this can take several hours. If available, load from file instead
# taxa_joined <- assignTaxonomy(seqtab_joined, unite.ref, multithread = MULTITHREAD, tryRC = TRUE)
# saveRDS(taxa_joined, "./data/sensitivity_taxa_joined.Rds")
taxa_joined <- readRDS("./data/sensitivity_taxa_joined.Rds")

# Combine elements into physeq
physeq_joined <- phyloseq(otu_table(seqtab_joined,taxa_are_rows=FALSE),
                          sample_data(sampledata),
                          tax_table(taxa_joined))

# Remove samples with zero total counts
physeq_joined_nonzero <- prune_samples(sample_sums(physeq_joined) > 0, physeq_joined)

# Ordinate
# ordination <- ordinate(physeq_joined_nonzero, "NMDS", "bray", k=2)
# saveRDS(ordination, "./data/sensitivity_ordination.Rds")
ordination <- readRDS("./data/sensitivity_ordination.Rds")
plot_ordination(physeq_joined_nonzero, ordination, type="samples",
                col="truncQ") +
  coord_cartesian(xlim=c(-0.096, -0.09), ylim=c(-0.002, 0.001)) +
  NULL

ordination_scores <- scores(ordination)
cor(ordination_scores[,1], as.matrix(sample_data(physeq_joined_nonzero)[,1:3]))
cor(ordination_scores[,2], as.matrix(sample_data(physeq_joined_nonzero)[,1:3]))
# No strong correlation with filtering parameters

which(ordination_scores[,1] > 15 & ordination_scores[,2] < 0)
which(ordination_scores[,1] > 15 & ordination_scores[,2] > 0)
# The outliers come from two samples: runB69PP_BMI_Plate4WellH3_ITS and runB69PP_BMI_Plate4WellH1_ITS
outliers <- c("runB69PP_BMI_Plate4WellH3_ITS","runB69PP_BMI_Plate4WellH1_ITS")

# Distances between samples within the same parameter set
hist(vegdist(seqtabs[[1]]))
hist(vegdist(seqtabs[[5]]))

# Distances between parameter sets for the same sample
par(mfrow=c(2,2))
hist(vegdist(seqtab_joined[grep("runB69PP_BMI_Plate13WellA12_ITS", rownames(seqtab_joined)),]))
hist(vegdist(seqtab_joined[grep("runC25G9_BMI_Plate70WellF7_ITS", rownames(seqtab_joined)),]))
hist(vegdist(seqtab_joined[grep("runC25G9_BMI_Plate70WellC10_ITS", rownames(seqtab_joined)),]))

# FINDING: Within the same parameter set, very little overlap between 
# samples (BCdist ~ 1). Across different parameter sets on the same sample,
# overlap can be large or small depending on the sample.

# Test this again, but aggregating to the species level, so that
# BCdist of samples with the same parameter sets is << 1.

physeq_joined_spp <- tax_glom(physeq_joined, taxrank="Species")
physeq_joined_spp_nonzero <- prune_samples(sample_sums(physeq_joined_spp) > 0, physeq_joined_spp)
# FYI the spp-aggregated sets don't contain the outlier samples, probably
# because those samples only contained taxa which were not classified to the
# spp-level.

physeqs_spp <- list()
for(i in 1:length(physeqs)) {
  physeqs_spp[[i]] <- tax_glom(physeqs[[i]], taxrank="Species")
}
names(physeqs_spp) <- param_sets

physeqs_spp_nonzero <- list()
for(i in 1:length(physeqs_spp)) {
  physeqs_spp_nonzero[[i]] <- prune_samples(sample_sums(physeqs_spp[[i]]) > 0, physeqs_spp[[i]])
}
names(physeqs_spp_nonzero) <- param_sets

# Plot ordination
set.seed(100101)
ordination_spp <- ordinate(physeq_joined_spp_nonzero, "NMDS", "bray", k=2)
plot(ordination_spp)

ordination_spp_scores <- scores(ordination_spp)
cor(ordination_spp_scores[,1], as.matrix(sample_data(physeq_joined_spp_nonzero)[,1:3]))
# No strong correlations with filtering parameters

which(ordination_spp_scores[,1] < -0.4)
# Outliers all come from runB69PP_BMI_Plate13WellA12_ITS

# Distances between samples in the same parameter set
lapply(physeqs_spp_nonzero, function(x) range(vegdist(otu_table(x)), na.rm = TRUE))
# as compared to non-agglomerated dataset
lapply(physeqs, function(x) range(vegdist(otu_table(x)), na.rm = TRUE))

# Distances between parameter sets for the same sample
vegdist(otu_table(physeq_joined_spp)[seq(from=1, to=125, by=5),])
hist(vegdist(otu_table(physeq_joined_spp)[seq(from=1, to=125, by=5),]))
# The finding re: distances is upheld, regardless of spp agglomeration

# Plot ordination of each sample across different parameter sets!
par(mfrow=c(2,3))
for(i in 1:5) {
  sample <- get_variable(physeq_joined_spp, "sample")[i]
  ps.subset <- subset_samples(physeq_joined_spp, sample=sample, sample_sums(physeq_joined_spp) > 0)
  ordination.subset <- ordinate(ps.subset, "NMDS", "bray", k=2)
  plot(ordination.subset, main=sample)
}


# Perform permANOVA
ps_joined_dist <- vegdist(otu_table(physeq_joined_nonzero))
adonis(ps_joined_dist ~ maxEE.F + truncQ,
       data = sample_data(physeq_joined_nonzero),
       strata = get_variable(physeq_joined_nonzero, "sample"),
       permutations=999)

# Try again with spp-agglomerated version
ps_joined_spp_dist <- vegdist(otu_table(physeq_joined_spp_nonzero))
adonis(ps_joined_spp_dist ~ maxEE.F + truncQ,
       data = sample_data(physeq_joined_spp_nonzero), 
       strata = get_variable(phseq_joined_spp_nonzero, "sample"),
       permutations=999)

# FINDING: Where there are compositional differences between parameter 
# sets, the differences are due to truncQ, not maxEE

