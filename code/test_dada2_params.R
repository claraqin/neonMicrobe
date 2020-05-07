# test_dada2_params.R
#
# DOWNLOAD WORKAROUND VERSION -- THIS TAKES AS INPUT A DIRECTORY
# CONTAINING FASTQ FILES TO BE TESTED UPON. YOU SHOULD STILL UPDATE
# params.R TO MATCH YOUR LOCAL FILE STRUCTURE, BUT YOU DO NOT NEED
# TO RUN 00_new_server_setup.R BEFORE RUNNING THIS SCRIPT. IN
# ADDITION, YOU WILL NEED TO DECOMPRESS THE 'C25G9_sample.zip' FILE,
# AVAILABLE AT THE LINK BELOW, AND PASTE ITS CONTENTS INTO 'PATH_CUT',
# DEFINED ON LINE ~57.
#
# LINK TO SAMPLE FILES: https://drive.google.com/file/d/1sIyp3_Ne1vF7qLf-racqoyNbEL_1qZGL/view?usp=sharing
#
# This script tests the effects of various filterAndTrim and dada2 parameters on
# - proportion / number of reads remaining in a sample
# - diversity estimates of a sample
# - taxonomic resolution of a sample


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

### TEMPORARY ############
# Create required directories
BASE_DIR <- PRESET_OUTDIR_SEQUENCE
if(!dir.exists(PRESET_OUTDIR_SEQUENCE)) dir.create(PRESET_OUTDIR_SEQUENCE, recursive=TRUE)
if(!dir.exists(PRESET_OUTDIR_SEQMETA)) dir.create(PRESET_OUTDIR_SEQMETA, recursive=TRUE)
if(!dir.exists(PRESET_OUTDIR_SOIL)) dir.create(PRESET_OUTDIR_SOIL, recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "ITS"))) dir.create(file.path(BASE_DIR, "ITS"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "16S"))) dir.create(file.path(BASE_DIR, "16S"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "ITS", "0_unzipped"))) dir.create(file.path(BASE_DIR, "ITS", "0_unzipped"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "ITS", "1_filtN"))) dir.create(file.path(BASE_DIR, "ITS", "1_filtN"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "ITS", "2_cutadapt"))) dir.create(file.path(BASE_DIR, "ITS", "2_cutadapt"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "ITS", "3_filtered"))) dir.create(file.path(BASE_DIR, "ITS", "3_unzipped"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "ITS", "4_seqtabs"))) dir.create(file.path(BASE_DIR, "ITS", "4_seqtabs"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "ITS", "track_reads"))) dir.create(file.path(BASE_DIR, "ITS", "track_reads"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "16S", "0_unzipped"))) dir.create(file.path(BASE_DIR, "16S", "0_unzipped"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "16S", "1_filtN"))) dir.create(file.path(BASE_DIR, "16S", "1_filtN"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "16S", "2_cutadapt"))) dir.create(file.path(BASE_DIR, "16S", "2_cutadapt"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "16S", "3_filtered"))) dir.create(file.path(BASE_DIR, "16S", "3_unzipped"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "16S", "4_seqtabs"))) dir.create(file.path(BASE_DIR, "16S", "4_seqtabs"), recursive=TRUE)
if(!dir.exists(file.path(BASE_DIR, "16S", "track_reads"))) dir.create(file.path(BASE_DIR, "16S", "track_reads"), recursive=TRUE)
### END TEMPORARY SECTION ########

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
  # regmatches(list.files(PATH_UNZIPPED), gregexpr("^run[A-Za-z0-9]*", list.files(PATH_UNZIPPED)))
  regmatches(list.files(PATH_CUT), gregexpr("^run[A-Za-z0-9]*", list.files(PATH_CUT)))
))

##############
# Track proportion of reads remaining after filterAndTrim under different params
i = 1
# This is an arbitrary selection; higher-quality runs include runBTJKN (#15), runC25G9 (#18), runC3CN4 (#19)
runID <- unique_runs[i]
runID
cutFs <- sort(list.files(PATH_CUT, pattern = paste0(runID, ".*_R1.fastq"), full.names = TRUE))
cutRs <- sort(list.files(PATH_CUT, pattern = paste0(runID, ".*_R2.fastq"), full.names = TRUE))

# To cut down on computation time, select 5 samples from the run:
if(length(cutFs) > 5) cutFs <- cutFs[round(quantile(1:length(cutFs)))]
if(length(cutRs) > 5) cutRs <- cutRs[round(quantile(1:length(cutRs)))]

# Plot quality profiles
for (sample_no in 1:length(cutFs)) {
  gridExtra::grid.arrange(plotQualityProfile(cutFs[[sample_no]]),
                          plotQualityProfile(cutRs[[sample_no]]), nrow=1)  
}

# Generalize previous version's setup to allow testing filterAndTrim on
# multiple sets of parameters at a time. The following tests 25 sets at a time:
grid <- expand.grid(2^seq(1,5), 2^seq(1,5))
params <- matrix(
  c(grid[,1],     # maxEE.F
    grid[,1],     # maxEE.R
    grid[,2],     # truncQ
    rep(50, nrow(grid))), # minLen (constant, for now)
  byrow=FALSE, ncol=4,
  dimnames = list(NULL, c("maxEE.F", "maxEE.R", "truncQ", "minLen"))
)
params

param_sets <- apply(params, 1, function(x) paste(c(rbind(c("maxEE.F", "maxEE.R", "truncQ", "minLen"), x)), collapse="_"))

# Create output filenames:
# Since there are 25 parameter sets, and 5 files in the sample,
# each list has 25 list elements, and each list element consists of 5 filenames.
filtFs <- lapply(param_sets, function(x) file.path(PATH_TEST, x, basename(cutFs)))
filtRs <- lapply(param_sets, function(x) file.path(PATH_TEST, x, basename(cutRs)))

# Run filterAndTrim on samples with all 25 parameter sets
out_list <- list()
for (i in 1:nrow(params)) {
  out <- filterAndTrim(
    fwd = cutFs, filt = filtFs[[i]], rev = cutRs, filt.rev = filtRs[[i]], 
    compress = TRUE, multithread = TRUE, maxN = 0,
    maxEE = params[i,1:2], truncQ = params[i,3], minLen = params[i,4]
  )
  out_list[[i]] <- as.data.frame(out) %>% mutate(prop.out = reads.out / reads.in)
}

names(out_list) <- param_sets

# Summarize and plot number/proportion of reads remaining
library(ggplot2)

# for disaggregated number of reads remaining
data.frame(
  maxEE = rep(params[,1], each=5),
  truncQ = rep(params[,3], each=5),
  reads.out = unlist(lapply(out_list, function(x) x[,"reads.out"]))
) %>%
  ggplot(aes(x=maxEE, y=truncQ, col=reads.out)) +
  geom_jitter() +
  scale_colour_gradient(low="blue", high="red")

# for median number of reads remaining
data.frame(
  maxEE = params[,1],
  truncQ = params[,3],
  reads.out = unlist(lapply(out_list, function(x) median(x[,"reads.out"])))
) %>%
  ggplot(aes(x=maxEE, y=truncQ, col=reads.out)) +
  geom_point() +
  scale_colour_gradient(low="blue", high="red")

# for median proportion of reads remaining
data.frame(
  maxEE = params[,1],
  truncQ = params[,3],
  prop.out = unlist(lapply(out_list, function(x) median(x[,"prop.out"])))
) %>%
  ggplot(aes(x=maxEE, y=truncQ, col=prop.out)) +
  geom_point() +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Median proportion of reads remaining\nafter filterAndTrim") +
  xlab("maxEE: after truncation, max allowable expected errors") +
  ylab("truncQ: quality score at which point to truncate read")
# Alternate version
data.frame(
  maxEE = params[,1],
  truncQ = params[,3],
  prop.out = unlist(lapply(out_list, function(x) median(x[,"prop.out"])))
) %>%
  ggplot(aes(x=truncQ, y=prop.out, col=maxEE)) +
  geom_point(alpha=0.6) +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Proportion of reads remaining\nafter filterAndTrim") +
  xlab("truncQ: quality score at which point to truncate read") +
  ylab("Prop. reads remaining (median of run)")

# These plots don't tell us much, because there is no guarantee that 
# any of these reads are of sufficient quality to be assigned taxonomy,
# or that they will merge.

###############
# Run parts of dada2 (post-filterAndTrim) to see downstream effects on merging, 
# taxonomic resolution, diversity

# Ensure that all of the output filenames exist
filtFs <- lapply(filtFs, function(x) { x <- x[file.exists(x)] })
filtRs <- lapply(filtRs, function(x) { x <- x[file.exists(x)] })

n_merged <- list()
prop_merged <- list()
seqtabs <- list()
dadaFs_list <- list()
dadaRs_list <- list()
system.time({
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
  
  ## TODO: Record intermediate metric here: the number of
  ##       sequence variants partitioned from the full set
  ##       of reads after the denoising algorithm.
  dadaFs_list[[i]] <- dadaFs
  dadaRs_list[[i]] <- dadaRs
  
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
  # hist(nchar(getSequences(seqtab.nochim)))
  
  seqtabs[[i]] <- seqtab.nochim
}
})
names(n_merged) <- param_sets
names(prop_merged) <- param_sets
names(dadaFs_list) <- param_sets
names(dadaRs_list) <- param_sets
names(seqtabs) <- param_sets

# Save the lists
saveRDS(dadaFs_list, "./data/sensitivity_dadaFs_list.Rds")
saveRDS(dadaRs_list, "./data/sensitivity_dadaRs_list.Rds")
saveRDS(seqtabs, "./data/sensitivity_seqtabs_list.Rds")

n_merged # This is the number of ASVs (as partitioned by dada) that successfully merged
prop_merged # This is proportion of ASVs (as partitioned by dada) that successfully merged
# The below is the proportion of forward reads that were assigned an ASV
prop_Fs_mapped_to_asv <- lapply(dadaFs_list, function(x) lapply(x, function(y) mean(!is.na(y$map))))
prop_Fs_mapped_to_asv_mat <- matrix(
  unlist(prop_Fs_mapped_to_asv),
  ncol=5,nrow=25,dimnames=list(param_sets, basename(cutFs))
)
head(prop_Fs_mapped_to_asv_mat)

# Rate of successfully merged ASVs per pre-filterAndTrim read
# (Not sure if this is meaningful -- treats denoising algorithm as a black box)
merged_variants_per_input_read <- lapply(n_merged, function(x) x/out_list[[1]][,"reads.in"])
merged_variants_per_input_read

# Summarize
summary_df <- 
  data.frame(
  maxEE = params[,1],
  truncQ = params[,3],
  prop.Fs.mapped = unlist(lapply(prop_Fs_mapped_to_asv, function(x) median(unlist(x)))),
  prop.merged = unlist(lapply(prop_merged, function(x) median(x))),
  rate.merged = unlist(lapply(merged_variants_per_input_read, function(x) median(x)))
)

# Plot the proportion of forward reads that were assigned to an ASV
summary_df %>%
  ggplot(aes(x=maxEE, y=truncQ, col=prop.Fs.mapped)) +
  geom_point() +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Proportion of (forward) reads assigned to an ASV") +
  xlab("maxEE: after truncation, max allowable expected errors") +
  ylab("truncQ: quality score at which point to truncate read")
# Alternate version
summary_df %>%
  ggplot(aes(x=truncQ, y=prop.Fs.mapped, col=maxEE)) +
  geom_point(alpha=0.6) +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Proportion of reads assigned to an ASV") +
  xlab("truncQ: quality score at which point to truncate read") +
  ylab("Prop. (forward) reads assigned to ASV\n(median of run)")

# Plot the proportion of ASVs that successfully merged
summary_df %>%
  ggplot(aes(x=maxEE, y=truncQ, col=prop.merged)) +
  geom_point() +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Proportion of ASVs that successfully merged") +
  xlab("maxEE: after truncation, max allowable expected errors") +
  ylab("truncQ: quality score at which point to truncate read")
# Alternate version
summary_df %>%
  ggplot(aes(x=truncQ, y=prop.merged, col=maxEE)) +
  geom_point(alpha=0.6) +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Proportion of ASVs that successfully merged") +
  xlab("truncQ: quality score at which point to truncate read") +
  ylab("Rate of successful merges (median of run)")

# Plot the rate of merged ASVs per pre-filter read
# (Not sure if this is meaningful -- treats denoising algorithm as a black box)
summary_df %>%
  ggplot(aes(x=maxEE, y=truncQ, col=rate.merged)) +
  geom_point() +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Successfully merged ASVs\nper input read") +
  xlab("maxEE: after truncation, max allowable expected errors") +
  ylab("truncQ: quality score at which point to truncate read")
# Alternate version
# summary_df %>%
#   ggplot(aes(x=maxEE, y=rate.merged)) +
#   geom_point() +
#   ggtitle("Successfully merged ASVs\nper input read") +
#   xlab("maxEE: after truncation, max allowable expected errors") +
#   ylab("Successfully merged ASVs per input read\n(median of run)")
summary_df %>%
  ggplot(aes(x=truncQ, y=rate.merged, col=maxEE)) +
  geom_point(alpha=0.6) +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Successfully merged ASVs\nper input read") +
  xlab("truncQ: quality score at which point to truncate read") +
  ylab("Successfully merged ASVs per input read\n(median of run)")

# These results indicate that the best parameters for ensuring
# a high number of successfully merged ASVs are low values
# for truncQ (2-8) and any values for maxEE (2-32)

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

saveRDS(taxas, "./data/sensitivity_taxas.Rds")


###################
# Evaluate effects on taxonomic resolution.
# What's the best way to evaluate taxonomic resolution?
# Number of taxa without a species-level assignment?
physeqs <- list()
abund_spp_classification <- list() # Total abundance of ASVs with species-level classifications
# prop_spp_classification <- list() # Proportion-by-abundance of ASVs with species-level classification
n_spp_classification <- list() # Number of unique ASVs with species-level classification
for(i in 1:length(seqtabs)) {
  physeqs[[i]] <- phyloseq(otu_table(seqtabs[[i]],taxa_are_rows=FALSE), tax_table(taxas[[i]]))
  abund_spp_classification[[i]] <- sum(otu_table(subset_taxa(physeqs[[i]], !is.na(Species))))
  # prop_spp_classification[[i]] <- abund_spp_classification[[i]] / sum(otu_table(physeqs[[i]]))
  n_spp_classification[[i]] <- sum(!is.na(tax_table(physeqs[[i]])[,"Species"]))
}
abund_spp_classification
# prop_spp_classification
n_spp_classification

# Remake summary_df with additional columns
summary_df <- 
  data.frame(
    maxEE = params[,1],
    truncQ = params[,3],
    prop.Fs.mapped = unlist(lapply(prop_Fs_mapped_to_asv, function(x) median(unlist(x)))),
    prop.merged = unlist(lapply(prop_merged, function(x) median(x))),
    rate.merged = unlist(lapply(merged_variants_per_input_read, function(x) median(x))),
    abund.spp.classification = unlist(abund_spp_classification),
    n.spp.classification = unlist(n_spp_classification)
  )

# Plot the abundance of species-level classifications in each physeq
# (Not sure if this is meaningful -- treats everything in between as black box)
summary_df %>%
  ggplot(aes(x=maxEE, y=truncQ, col=abund.spp.classification)) +
  geom_point() +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Abundance of ASVs with species-level\nclassifications") +
  xlab("maxEE: after truncation, max allowable expected errors") +
  ylab("truncQ: quality score at which point to truncate read")

# Two-dim versions
summary_df %>%
  ggplot(aes(x=maxEE, y=abund.spp.classification)) +
  geom_point(alpha=0.6) +
  ggtitle("Abundance of ASVs with species-level\nclassifications") +
  xlab("maxEE: after truncation, max allowable expected errors") +
  ylab("Total abundance of spp-level ASVs")
summary_df %>%
  ggplot(aes(x=truncQ, y=abund.spp.classification, col=maxEE)) +
  geom_point(alpha=0.6) +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Abundance of ASVs with species-level\nclassifications") +
  xlab("truncQ: quality score at which point to truncate read") +
  ylab("Total abundance of spp-level ASVs")

# Plot the number of unique ASVs with species-level classification
summary_df %>%
  ggplot(aes(x=maxEE, y=truncQ, col=n.spp.classification)) +
  geom_point() +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Number of unique ASVs with\nspecies-level classifications") +
  xlab("maxEE: after truncation, max allowable expected errors") +
  ylab("truncQ: quality score at which point to truncate read")

# Two-dim versions
summary_df %>%
  ggplot(aes(x=maxEE, y=n.spp.classification)) +
  geom_point(alpha=0.6) +
  ggtitle("Number of unique ASVs with\nspecies-level classifications") +
  xlab("maxEE: after truncation, max allowable expected errors") +
  ylab("Total abundance of spp-level ASVs")
summary_df %>%
  ggplot(aes(x=truncQ, y=n.spp.classification, col=maxEE)) +
  geom_point(alpha=0.6) +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Number of unique ASVs with\nspecies-level classifications") +
  xlab("truncQ: quality score at which point to truncate read") +
  ylab("No. unique spp-level ASVs")

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
estimate_richness(physeqs[[1]], measures=c("Shannon"))
estimate_richness(physeqs[[1]], measures=c("Shannon"), split=FALSE) # aggregate into one measure

shannon_list <- list()
for(i in 1:length(physeqs)) {
  shannon_list[[i]] <- suppressWarnings(estimate_richness(physeqs[[i]], measures=c("Shannon"), split=FALSE)[1,1])
}

# Remake summary_df with additional columns (Shannon index)
summary_df <- 
  data.frame(
    maxEE = params[,1],
    truncQ = params[,3],
    prop.Fs.mapped = unlist(lapply(prop_Fs_mapped_to_asv, function(x) median(unlist(x)))),
    prop.merged = unlist(lapply(prop_merged, function(x) median(x))),
    rate.merged = unlist(lapply(merged_variants_per_input_read, function(x) median(x))),
    abund.spp.classification = unlist(abund_spp_classification),
    n.spp.classification = unlist(n_spp_classification),
    shannon.div = unlist(shannon_list)
  )
summary_df %>%
  ggplot(aes(x=truncQ, y=shannon.div, col=maxEE)) +
  geom_point(alpha=0.6) +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Shannon diversity of aggregated samples") +
  xlab("truncQ: quality score at which point to truncate read") +
  ylab("Shannon diversity (H)")
