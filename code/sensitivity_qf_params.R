# sensitivity_qf_params.R
#
# This script tests the effects of various filterAndTrim and dada2 parameters on
# - proportion / number of reads remaining in a sample
# - diversity estimates of a sample
# - taxonomic resolution of a sample
# ... for 16S sequences
#
# This script was based on test_dada2_params.R, but was revised to process 16S sequences,
# and to test different parameters.
#
# Parameters to test: maxEE of R2 reads and truncLen of R2 reads.

# TODO: Make the testable parameters more flexible. I (Clara) attempted to do this at first
#       by replacing "maxEE.R" and "truncLen" with PARAM1 and PARAM2 respectively, but this
#       isn't easily passed to functions that use non-standard evaluation, like ggplot aes.
#       (This might be addressed by the new tidy evaluation semantics of dplyr, which makes it
#       so we can use enquo() and !!.)
#       Flexibility is also fundamentally limited in that the user of this code must know
#       where in the processing pipeline to input the varying parameters.

##############
# Preset values for sensitivity analysis

# Quality filter parameters that we do not vary in this sensitivity analysis
MAX_EE_FWD <- 8
TRUNC_LEN_FWD <- 240
MIN_LEN <- 50

# Parameters specific to this sensitivity analysis
N_SAMPLES <- 40
DIRNAME_TEST <- "qf_test_11-09-2020"
runIDs <- c("C5B2R", "CBTWG")

# Parameter value grid. The following allows testing of two quality filtering parameters at a time.
PARAM1 <- "maxEE.R"
PARAM2 <- "truncLen.R"
grid <- expand.grid(c(4, 8, 16), c(160, 220, 280))
params <- matrix(
  c(grid[,1],     # PARAM1
    grid[,2]),    # PARAM2
  byrow=FALSE, ncol=2,
  dimnames = list(NULL, c(PARAM1, PARAM2))
)
param_sets <- apply(params, 1, function(x) paste(c(rbind(c(PARAM1, PARAM2), x)), collapse="_"))


##############
# Load parameters, libraries, and filepaths

# Load parameters from params.R, tools from utils.R
source("./code/params.R")
source("./R/utils.R")

# Load libraries
library(dada2)
library(ShortRead)
library(Biostrings)
library(tibble)
library(dplyr)
library(vegan)
library(phyloseq)
library(ggplot2)
library(tidyr)

# Generate filepath names
if(is.null(PRESET_OUTDIR_SEQUENCE) | PRESET_OUTDIR_SEQUENCE == "") {
  PATH_16S <- file.path(PRESET_OUTDIR, "raw_sequence", "16S")
} else {
  PATH_16S <- file.path(PRESET_OUTDIR, PRESET_OUTDIR_SEQUENCE, "16S")
}
PATH_RAW <- file.path(PATH_16S, "0_raw")
PATH_TRIMMED <- file.path(PATH_16S, "1_trimmed")
PATH_TEST <- file.path(PATH_16S, DIRNAME_TEST)
dir.create(file.path(PATH_TEST, "results"))

##############
# Runs to use:
rawFs <- sort(list.files(PATH_RAW, pattern = paste0("(",paste0("(", runIDs, ")", collapse="|"),")", ".*_R1\\.fastq"), full.names = TRUE))
rawRs <- sort(list.files(PATH_RAW, pattern = paste0("(",paste0("(", runIDs, ")", collapse="|"),")", ".*_R2\\.fastq"), full.names = TRUE))
# cutFs <- sort(list.files(PATH_TRIMMED, pattern = paste0("(",paste0("(", runIDs, ")", collapse="|"),")", ".*_R1\\.fastq"), full.names = TRUE))
# cutRs <- sort(list.files(PATH_TRIMMED, pattern = paste0("(",paste0("(", runIDs, ")", collapse="|"),")", ".*_R2\\.fastq"), full.names = TRUE))

# Remove any files that only have forward or reverse reads
matched_fn <- remove_unmatched_files(rawFs, rawRs)
rawFs <- matched_fn[[1]]
rawRs <- matched_fn[[2]]

# To cut down on computation time, select up to N_SAMPLES samples from the runs, up to N_SAMPLES/length(runIDs) from each run.
if(length(rawFs) > N_SAMPLES) {
  rawFs_subset <- c()
  rawRs_subset <- c()
  for(i in 1:length(runIDs)) {
    rawFs_runID <- rawFs[grep(runIDs[i], rawFs)]
    rawRs_runID <- rawRs[grep(runIDs[i], rawRs)]
    if(length(rawFs_runID) > 5) {
      rawFs_subset <- c(rawFs_subset, rawFs_runID[round(quantile(1:length(rawFs_runID), probs=seq(0,1,length.out=N_SAMPLES/length(runIDs))))])
      rawRs_subset <- c(rawRs_subset, rawRs_runID[round(quantile(1:length(rawRs_runID), probs=seq(0,1,length.out=N_SAMPLES/length(runIDs))))])
    } else {
      rawFs_subset <- c(rawFs_subset, rawFs_runID)
      rawRs_subset <- c(rawRs_subset, rawRs_runID)
    }
  }
  rawFs <- rawFs_subset
  rawRs <- rawRs_subset
}

# Plot quality profiles
gridExtra::grid.arrange(plotQualityProfile(rawFs, aggregate=TRUE),
                        plotQualityProfile(rawRs, aggregate=TRUE))

# Split complete filenames into "basenames" and directory names
fn_base <- basename(c(rawFs, rawRs))
PATH_PARAMSETS <- file.path(PATH_TEST, param_sets)

t1 <- Sys.time()
# Trim reads based on the primer sequences supplied in params.r
trim_trackReads <- trimPrimers16S(fn_base, PATH_RAW, PATH_TRIMMED, "CCTACGGGNBGCASCAG", "GACTACNVGGGTATCTAATCC")

# Run quality filter on 16S sequences
t2 <- Sys.time()
filter_trackReads <- list()
for(i in 1:length(param_sets)) {
  filter_trackReads[[i]] <- qualityFilter16S(
    fn_base,
    PATH_TRIMMED,
    PATH_PARAMSETS[[i]],
    maxEE=c(MAX_EE_FWD, params[i,1]), # Vary maxEE.R
    truncLen=c(TRUNC_LEN_FWD, params[i,2]), # Vary truncLen.R
    minLen=MIN_LEN,
    multithread=MULTITHREAD)
  rownames(filter_trackReads[[i]]) <- paste0(param_sets[i], "|", rownames(filter_trackReads[[i]]))
}
filter_trackReads_mat <- do.call(rbind, filter_trackReads)

###############
# Run parts of dada2 (post-filterAndTrim) to see downstream effects on merging,
# taxonomic resolution, diversity

seqtabs <- dada_trackReads <- lapply(1:length(param_sets), function(x) lapply(1:length(runIDs), function(y) list()))

t3 <- Sys.time() # t1 - t3 on 11/09/2020: 8 min

# system.time({ # When run on runC5B2R_BMI_Plate77WellA10_16S_C5B2R and runC5B2R_BMI_Plate77WellA11_16S_C5B2R,
#               # on 8 cores, takes ~35 min.
for(i in 1:length(param_sets)) {
  for(j in 1:length(runIDs)) {
    message("Sensitivity analysis: parameter set ", param_sets[i], ", sequencing run ",  runIDs[j])

    # Retrieve only those files associated with the appropriate parameter set and runID
    fn_base.star <- grep(runIDs[j], fn_base, value=TRUE)

    seqtab.list <- runDada16S(fn_base.star, PATH_PARAMSETS[i], MULTITHREAD, VERBOSE, seed=11001100)
    seqtabs[[i]][[j]] <- seqtab.list$seqtab.nochim
    dada_trackReads[[i]][[j]] <- seqtab.list$track

    rownames(seqtabs[[i]][[j]]) <- paste0(param_sets[i], "|", rownames(seqtabs[[i]][[j]]))
    rownames(dada_trackReads[[i]][[j]]) <- paste0(param_sets[i], "|", rownames(dada_trackReads[[i]][[j]]))
  }
}
# })
dada_trackReads_mat <- do.call(rbind, lapply(dada_trackReads, function(x) do.call(rbind, x)))

t4 <- Sys.time()
# Combine all tracking tables -- first requires modification of trim_trackReads
trim_trackReads_mat <- do.call(rbind, lapply(1:length(param_sets), function(i) {
  rownames(trim_trackReads) <- paste0(param_sets[i], "|", rownames(trim_trackReads))
  trim_trackReads
}))

track <- Reduce(
  function(x, y, ...) transform(merge(x, y, by=0, all = TRUE, ...), row.names=Row.names, Row.names=NULL),
  list(trim_trackReads_mat[,1],
       filter_trackReads_mat,
       dada_trackReads_mat)
)
names(track) <- c("input", "trimmed", "filtered", colnames(dada_trackReads_mat))
track[is.na(track)] <- 0

# For community analyses, join sequencing runs together for each parameter set
seqtabs_joinrun <- lapply(1:length(param_sets), function(x) list())
for(i in 1:length(param_sets)) {
  seqtabs_joinrun[[i]] <- mergeSequenceTables(tables = seqtabs[[i]])
}

# Save the data so far
saveRDS(seqtabs_joinrun, file.path(PATH_TEST, "results", "sensitivity_seqtabs_joinrun_list.Rds"))
write.csv(track, file.path(PATH_TEST, "results", "sensitivity_trackReads.csv"), row.names=TRUE)

# seqtabs_joinrun <- readRDS(file.path(PATH_TEST, "results", "sensitivity_seqtabs_joinrun_list.Rds"))

####
# Plot read counts at each step

# Parse parameter values from the rownames of a data frame
track <- parseParamsFromRownames(track, PARAM1, PARAM2)

# Reshape read tracking table
track_long <- tidyr::gather(track, key = "step", value = "reads", input:nonchim)
# Exclude metrics associated only with forward reads
track_long <- track_long[!grepl("F$", track_long$step),]
# Aggregate read counts by run ID
track_long[["step"]] <- factor(track_long[["step"]], levels=colnames(track)[1:9])
track_aggRun <- group_by(track_long, maxEE.R, truncLen.R, runID, step) %>%
  dplyr::summarise(reads = sum(reads))

# Plot
theme_set(theme_bw())
ggplot(track_aggRun, aes(x=step, y=reads, col=as.factor(maxEE.R))) +
  geom_line(aes(linetype=as.factor(truncLen.R), group=interaction(maxEE.R, truncLen.R))) +
  facet_wrap(~runID, scales="free_y") +
  labs(linetype="truncLen.R", color="maxEE.R") +
  scale_linetype_manual(values=c("dotted", "dashed", "solid")) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(file.path(PATH_TEST, "results", "track_reads_plot.png"), width=9, height=3.5, units="in")

###
# Assign taxonomy using the UNITE database
# Warning: this can take several days. If available, load from file
t5 <- Sys.time()
taxas <- list()
tax_t0 <- Sys.time()
tax_t <- c()
for (i in 1:length(seqtabs_joinrun)) { # When run on runC5B2R_BMI_Plate77WellA10_16S_C5B2R and runC5B2R_BMI_Plate77WellA11_16S_C5B2R,
  # on 8 cores, takes ~40 min.
  taxas[[i]] <- assignTaxonomy(seqtabs_joinrun[[i]], UNITE_REF_PATH, multithread = MULTITHREAD, tryRC = TRUE)
  tax_t <- c(tax_t, Sys.time())
}
names(taxas) <- param_sets
t6 <- Sys.time()

# taxas.print <- taxas[[1]]
# rownames(taxas.print) <- NULL
# tail(taxas.print)

# Save the data thus far
saveRDS(taxas, file.path(PATH_TEST, "results", "sensitivity_taxas_list.Rds"))
# taxas <- readRDS(file.path(PATH_TEST, "results", "sensitivity_taxas_list.Rds"))


###################
# Evaluate effects on taxonomic resolution.
# What's the best way to evaluate taxonomic resolution?
# Number of taxa without a species-level assignment?

physeqs <- list() # Create a list of physeq objects
abund_spp_classification <- list() # Total abundance of ASVs with species-level classifications
prop_spp_classification <- list() # Proportion-by-abundance of ASVs with species-level classification
n_spp_classification <- list() # Number of unique ASVs with species-level classification
for(i in 1:length(param_sets)) {
  physeqs[[i]] <- phyloseq(otu_table(seqtabs_joinrun[[i]],taxa_are_rows=FALSE), tax_table(taxas[[i]]))
  abund_spp_classification[[i]] <- tryCatch({
    sum(otu_table(subset_taxa(physeqs[[i]], !is.na(Species))))
  }, error = function(e) 0)
  prop_spp_classification[[i]] <- abund_spp_classification[[i]] / sum(otu_table(physeqs[[i]]))
  n_spp_classification[[i]] <- sum(!is.na(tax_table(physeqs[[i]])[,"Species"]))
}

# Total number of reads, regardless of assignment to taxonomy
n_reads <- lapply(physeqs, function(ps) sum(otu_table(ps)))
# Number of reads assigned to taxonomic rank. A two-tiered list:
# - First tier: minimum taxonomic rank. [1:7]
# - Second tier: parameter set. [1:9]
tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
n_reads_assigned_rank <- lapply(
  tax_ranks,
  function(rank) lapply(physeqs, function(ps) {
      tax <- data.frame(as(tax_table(ps), "matrix"))
      ind <- !is.na(tax[[rank]])
      if(all(ind==FALSE)) return(0)
      sum(otu_table(ps)[, ind])
    }
  )
)
n_reads_assigned_rank_df <- as.data.frame(do.call(cbind, lapply(n_reads_assigned_rank, unlist)))
# p_reads_assigned_rank_df <- as.data.frame(do.call(cbind, lapply(n_reads_assigned_rank, unlist)) / unlist(n_reads))
rownames(n_reads_assigned_rank_df) <- param_sets
colnames(n_reads_assigned_rank_df) <- tax_ranks
n_reads_assigned_rank_df <- parseParamsFromRownames(n_reads_assigned_rank_df, PARAM1, PARAM2) %>%
  dplyr::select(maxEE.R, truncLen.R, Kingdom:Species)
n_reads_assigned_rank_df

# Reshape assigned reads table, and calculate proportion
tidyr::gather(n_reads_assigned_rank_df, key="rank", value="n_assigned_reads", Kingdom:Species) %>%
  mutate(p_assigned_reads = n_assigned_reads / unlist(n_reads),
         rank = factor(rank, levels=tax_ranks)) ->
  reads_assigned_rank_df

# Plot
theme_set(theme_bw())
ggplot(reads_assigned_rank_df, aes(x=rank, y=p_assigned_reads, col=as.factor(maxEE.R))) +
  geom_line(aes(linetype=as.factor(truncLen.R), group=interaction(maxEE.R, truncLen.R))) +
  labs(linetype="truncLen.R", color="maxEE.R") +
  scale_linetype_manual(values=c("dotted", "dashed", "solid")) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(file.path(PATH_TEST, "results", "reads_assigned_taxonomy.png"), width=5, height=3.5, units="in")


##################
# Effects on alpha-diversity

# What's the best way to evaluate effects on diversity estimates?
# Do the parameters affect the shape of the rarefaction curves (evenness),
# or just the sequencing depth (number of reads)?
# # par(mfrow=c(2,2))
# # lapply(seqtabs_joinrun, function(x) {
# #   rarecurve(x, step=50, label=FALSE)
# # })
# # par(mfrow=c(1,1))
# # rarecurve(seqtabs[[4]], sample=50)
#
# # rarefactions <- list()
# # for(i in 1:length(param_sets)) {
# #   max_depth <- max(rowSums(seqtabs_joinrun[[i]]))
# #   rarefy(seqtabs_joinrun[[i]], sample=)
# # }
# # rarefy(seqtabs_joinrun[[1]], sample=c(50, 100, 150, 200))

t7 <- Sys.time()

# Use phyloseq's estimate_richness function
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

names(obsrich_list) <- names(shannon_list) <- param_sets
as.data.frame(cbind(obsrich_list, shannon_list)) %>%
  dplyr::rename(`Observed richness` = "obsrich_list",
                `Shannon index` = "shannon_list") %>%
  parseParamsFromRownames(PARAM1, PARAM2) %>%
  tidyr::gather(key="index", value="value", `Observed richness`:`Shannon index`) %>%
  mutate(value = unlist(value)) ->
  diversity_df


theme_set(theme_bw())
ggplot(diversity_df, aes(y=value, col=as.factor(maxEE.R), x=as.factor(truncLen.R))) +
  geom_boxplot() +
  facet_grid(index~., scales="free_y") +
  labs(col="maxEE.R", x="truncLen.R")
ggsave(file.path(PATH_TEST, "results", "alpha_diversity_plot.png"), width=5, height=5, units="in")


##############
# Evaluate effects on beta-diversity inference,
# using ordinations and permANOVA

t8 <- Sys.time()

# Combine all seqtabs into ONE physeq
seqtab_joined <- mergeSequenceTables(tables=seqtabs_joinrun)

# Sample data (parameters)
sampledata <- parseParamsFromRownames(seqtab_joined, PARAM1, PARAM2, keep_rownames=TRUE, keep_original_cols = FALSE)

# Taxa table
# Warning: this can take several hours. If available, load from file instead
t9 <- Sys.time()
# system.time({ # When run on runC5B2R_BMI_Plate77WellA10_16S_C5B2R and runC5B2R_BMI_Plate77WellA11_16S_C5B2R,
  # with 9 parameter sets, on 8 cores, takes ~15 min.
taxa_joined <- assignTaxonomy(seqtab_joined, UNITE_REF_PATH, multithread = MULTITHREAD, tryRC = TRUE)
# })

saveRDS(taxa_joined, file.path(PATH_TEST, "results", "sensitivity_taxa_joined.Rds"))
# taxa_joined <- readRDS("./data/sensitivity_taxa_joined.Rds")

t10 <- Sys.time()

# Combine elements into physeq
physeq_joined <- phyloseq(otu_table(seqtab_joined,taxa_are_rows=FALSE),
                          sample_data(sampledata),
                          tax_table(taxa_joined))

# Remove samples with zero total counts
physeq_joined_nonzero <- prune_samples(sample_sums(physeq_joined) > 0, physeq_joined)

# Optional: Keep only taxa with at least 10 reads
# physeq_joined_nonzero <- prune_taxa(taxa_sums(physeq_joined) > 10, physeq_joined)

# Ordinate
ordination <- ordinate(physeq_joined_nonzero, "NMDS", "bray", k=2)

saveRDS(ordination, file.path(PATH_TEST, "results", "sensitivity_ordination.Rds"))

t11 <- Sys.time()

ordination <- readRDS(file.path(PATH_TEST, "results", "sensitivity_ordination.Rds"))

theme_set(theme_bw())
sample_data(physeq_joined_nonzero)[,"maxEE.R.factor"] <- as.factor(get_variable(physeq_joined_nonzero, "maxEE.R"))
plot_ordination(physeq_joined_nonzero, ordination, type="samples",
                col="truncLen.R", shape="maxEE.R.factor") +
  # coord_cartesian(xlim=c(-0.096, -0.09), ylim=c(-0.002, 0.001)) +
  NULL
ggsave(file.path(PATH_TEST, "results", "sensitivity_ordination.png"), width=5, height=3.5, units="in")

ordination_scores <- scores(ordination)




####### PAUSED HERE ON 12/14/2020


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

