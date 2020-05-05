# test_dada2_params.R
#
# DOWNLOAD WORKAROUND VERSION -- THIS TAKES AS INPUT A DIRECTORY
# CONTAINING FASTQ FILES TO BE TESTED UPON. YOU SHOULD STILL UPDATE
# params.R TO MATCH YOUR LOCAL FILE STRUCTURE, BUT YOU DO NOT NEED
# TO RUN 00_new_server_setup.R BEFORE RUNNING THIS SCRIPT. IN
# ADDITION, YOU WILL NEED TO DECOMPRESS THE 'C25G9_sample.zip' FILE,
# AVAILABLE AT THE LINK BELOW, AND PASTE ITS CONTENTS INTO 'PATH_CUT',
# DEFINED ON LINE 61.
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
list_out <- list()
for (i in 1:nrow(params)) {
  out <- filterAndTrim(
    fwd = cutFs, filt = filtFs[[i]], rev = cutRs, filt.rev = filtRs[[i]], 
    compress = TRUE, multithread = TRUE, maxN = 0,
    maxEE = params[i,1:2], truncQ = params[i,3], minLen = params[i,4]
  )
  list_out[[i]] <- as.data.frame(out) %>% mutate(prop.out = reads.out / reads.in)
}

names(list_out) <- param_sets

# Summarize and plot number/proportion of reads remaining
library(ggplot2)

# for disaggregated number of reads remaining
data.frame(
  maxEE = rep(params[,1], each=5),
  truncQ = rep(params[,3], each=5),
  reads.out = unlist(lapply(list_out, function(x) x[,"reads.out"]))
) %>%
  ggplot(aes(x=maxEE, y=truncQ, col=reads.out)) +
  geom_jitter() +
  scale_colour_gradient(low="blue", high="red")

# for median number of reads remaining
data.frame(
  maxEE = params[,1],
  truncQ = params[,3],
  reads.out = unlist(lapply(list_out, function(x) median(x[,"reads.out"])))
) %>%
  ggplot(aes(x=maxEE, y=truncQ, col=reads.out)) +
  geom_point() +
  scale_colour_gradient(low="blue", high="red")

# for median proportion of reads remaining
data.frame(
  maxEE = params[,1],
  truncQ = params[,3],
  prop.out = unlist(lapply(list_out, function(x) median(x[,"prop.out"])))
) %>%
  ggplot(aes(x=maxEE, y=truncQ, col=prop.out)) +
  geom_point() +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Median proportion of reads remaining\nafter filterAndTrim") +
  xlab("maxEE: after truncation, max allowable expected errors") +
  ylab("truncQ: quality score at which point to truncate read")

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
names(n_merged) <- param_sets
names(prop_merged) <- param_sets
names(seqtabs) <- param_sets

n_merged # This is the number of ASVs (as partitioned by dada) that successfully merged
prop_merged # This is proportion of ASVs (as partitioned by dada) that successfully merged

# Rate of successfully merged ASVs per pre-filterAndTrim read
# (Not sure if this is meaningful -- treats denoising algorithm as a black box)
merged_variants_per_input_read <- lapply(n_merged, function(x) x/list_out[[1]][,"reads.in"])
merged_variants_per_input_read

# Summarize and plot the proportion of ASVs that successfully merged
data.frame(
  maxEE = params[,1],
  truncQ = params[,3],
  prop.merged = unlist(lapply(prop_merged, function(x) median(x)))
) %>%
  ggplot(aes(x=maxEE, y=truncQ, col=prop.merged)) +
  geom_point() +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Median proportion of ASVs successfully merged") +
  xlab("maxEE: after truncation, max allowable expected errors") +
  ylab("truncQ: quality score at which point to truncate read")

# Summarize and plot the rate of merged ASVs per pre-filter read
# (Not sure if this is meaningful -- treats denoising algorithm as a black box)
data.frame(
  maxEE = params[,1],
  truncQ = params[,3],
  rate.merged = unlist(lapply(merged_variants_per_input_read, function(x) median(x)))
) %>%
  ggplot(aes(x=maxEE, y=truncQ, col=rate.merged)) +
  geom_point() +
  scale_colour_gradient(low="blue", high="red") +
  ggtitle("Median rate of successfully merged ASVs\nper input read") +
  xlab("maxEE: after truncation, max allowable expected errors") +
  ylab("truncQ: quality score at which point to truncate read")

# These results indicate that the best parameters for ensuring
# a high number of successfully merged ASVs are low values
# for truncQ (2-8) and any values for maxEE (2-32)

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

