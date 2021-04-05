# neonMicrobe

`neonMicrobe` is a suite of functions for downloading, pre-processing, and assembling heterogeneous data around the NEON soil microbe marker gene sequence data. To do so, `neonMicrobe` downloads NEON data products from the NEON Data API and processes sequences using the DADA2 workflow. In addition, `neonMicrobe` offers a processing-batch infrastructure to encourage explicit versioning of processed data.

## Installation

The development version of `neonMicrobe` can be installed directly from this GitHub repo using this code:

```
install.packages("devtools")
devtools::install_github("claraqin/neonMicrobe")
```

## User-installed dependencies

In addition to the R package dependencies which are installed alongside `neonMicrobe`, users may also need to complete the following requirements before using some functions in `neonMicrobe`:

1. **For taxonomic assignment via DADA2**, you will need to install the latest taxonomic reference datasets for ITS or 16S sequences. Consult the [DADA2 taxonomic reference data webpage](https://benjjneb.github.io/dada2/training.html) for more information. For organizational purposes, we recommend keeping these files in the `data/tax_ref` subdirectory that is created after you run `makeDataDirectories()` (see "Input data" below).
2. **For trimming of ITS sequence primers**, you will need to install `cutadapt`. Installation instructions can be found [here](https://cutadapt.readthedocs.io/en/stable/installation.html). Once installed, you can tell `neonMicrobe` where to look for it by specifying the `cutadapt_path` argument each time you use the `trimPrimerITS` function. For an example, see the "Process 16S Sequences" vignette or the "Process ITS Sequences" vignette.

## Overview

Tutorials for `neonMicrobe` are available in the `vignettes` directory.

1. **Download NEON Data** – how to use the functions in this package to download the specific scope NEON soil microbe marker gene sequence data and associated data relevant to your analysis. Leverages the `neonUtilities` [R package](https://github.com/NEONScience/NEON-utilities).
2. **Process Sequences** – how to use the functions in this package to add associated environmental variables to the ASV tables. Leverages the `dada2` [R package](https://github.com/benjjneb/dada2). The dada denoising algorithm partitions reads into amplicon sequence variants (ASVs), which are finer in resolution than OTUs.
3. **Add Environmental Variables** – how to use the functions in this package to add associated environmental variables to the ASV tables. Joins the data together in the form of one or more Phyloseq objects.
4. (Optional) **Sensitivity Analysis** – how to use the functions in this package to test the effects of quality filtering parameters and decisions on the resulting ecological inference. Can be used to test for an acceptable range of custom parameters.

![NEON Ecosphere MS Figure-Making Workspace (3)](https://user-images.githubusercontent.com/12421420/111393342-ce937d00-8675-11eb-8b63-530aced18352.png)

## Data storage structure

### Input data

The **Download NEON Data** vignette demonstrates how to download NEON data, optionally writing to the file system. By default, the input data is downloaded into the following structure, which is created in the working directory after running `makeDataDirectories()`:

![NEON Ecosphere MS Figure-Making Workspace (10)](https://user-images.githubusercontent.com/12421420/113089173-f3badc00-919b-11eb-84e6-b7f9a2abbb72.png)

The tree structure in the upper-left represents the data directory structure constructed within the project root directory. Red dotted lines represent explicit linkages between NEON data products via shared data fields. (a) Sequence metadata is downloaded from NEON data product DP1.10108.001 (Soil microbe marker gene sequences) using the downloadSequenceMetadata() function. (b) Raw microbe marker gene sequence data is downloaded from NEON based on the sequence metadata using the downloadRawSequenceData() function. (c) Soil physical and chemical data is downloaded from NEON data product DP1.10086.001 using the downloadRawSoilData() function. (d) Taxonomic reference datasets (e.g. SILVA, UNITE) are added separately by the user.

### Output data

The **Process Sequences** and **Add Environmental Variables** vignettes demonstrate how to process the NEON data inputs into useful sample-abundance tables with accompanying environmental data.

(Under construction.)


## Methods paper

A methods paper describing the use of `neonMicrobe` is currently in review.

# Quick Start

```
# analyze-neon-greatplains-16s.R

# Download and analyze NEON soil microbe marker gene sequence data from select grassland sites
# for demonstration purposes

library(plyr)
library(dplyr)
library(neonUtilities)
library(neonMicrobe)

# On load, neonMicrobe says:
# neonMicrobe relies on a pre-generated directory structure. If this 
# is your first time using neonMicrobe, or you want to create a new 
# directory structure, first set your working directory to the location 
# where you would like to create this structure, and then run 
# makeDataDirectories(). Learn more in the vignette 'Download NEON Data'.

# neonMicrobe generates default output directories based on the current 
# working directory. To hold this constant, set a 'base directory' using
# setBaseDirectory().

setBaseDirectory()
makeDataDirectories()

meta_16s <- downloadSequenceMetadata(startYrMo = "2017-07", endYrMo = "2017-07",
                                     sites = c("KONZ", "CPER", "NOGP"),
                                     targetGene = "16S")

meta_16s_qc <- qcMetadata(meta_16s, pairedReads = "Y")

with(meta_16s_qc, table(siteID, sequencerRunID))

downloadRawSequenceData(meta_16s_qc)


# Process sequences
# See the Process ITS Sequences vignette for a more detailed explanation

library(ShortRead)
library(Biostrings)
library(dada2)

meta <- meta_16s_qc

fl_nm <- meta$rawDataFileName
head(fl_nm)

trim_trackReads <- trimPrimers16S(
  fl_nm, in_subdir = "raw", out_subdir = "1_trimmed",
  meta = meta
)

filter_trackReads <- qualityFilter16S(
  fl_nm, in_subdir = "1_trimmed", out_subdir = "2_filtered",
  meta = meta, truncLen = 220, maxEE = 8, multithread = TRUE
)

unique_runs <- unique(meta$sequencerRunID)

for(i in 1:length(unique_runs)) {
  meta_thisrun <- meta[which(meta$sequencerRunID==unique_runs[i]),]
  fl_nm_thisrun <- meta_thisrun$rawDataFileName
  dada_out <- runDada16S(
    fl_nm_thisrun, in_subdir = "2_filtered", meta = meta,
    out_seqtab = file.path(NEONMICROBE_DIR_OUTPUTS(), "greatplains",
                           paste0("greatplains_asv_", unique_runs[i], ".Rds")),
    out_track = file.path(NEONMICROBE_DIR_OUTPUTS(), "greatplains",
                          paste0("greatplains_track_", unique_runs[i], ".csv")),
    verbose = FALSE,
    multithread = TRUE
  )
}

seqtab_joined <- mergeSequenceTables(
  tables = file.path(NEONMICROBE_DIR_OUTPUTS(), "greatplains",
                     paste0("greatplains_asv_", unique_runs, ".csv"))
)


# Collapses NEON 16S ASV table

seqtab_collapse_filename <- file.path(NEONMICROBE_DIR_OUTPUTS(), "greatplains",
                                      "NEON_16S_seqtab_nochim_grasslands_COLLAPSED.Rds")
t0 <- Sys.time()
seqtab_collapse <- collapseNoMismatch(seqtab_joined)
saveRDS(seqtab_collapse, seqtab_collapse_filename)
t1 <- Sys.time()
t1 - t0 # Took 6.93 hours on socs-stats.ucsc.edu

seqtab_collapse <- readRDS(seqtab_collapse_filename)

# Add soil data

soils <- downloadSoilData(startYrMo = "2017-07", endYrMo = "2017-07",
                          sites = c("KONZ", "CPER", "NOGP"))

sampledata <- meta_16s_qc %>% distinct(dnaSampleID, .keep_all = T) %>%
  dplyr::select(sequencerRunID, dnaSampleID, dataQF.rawFiles, internalLabID.seq,
         siteID, collectDate, plotID, deprecatedVialID, geneticSampleID)
sampledata$date_ymd <- as.Date(format(as.Date(sampledata$collectDate, format="%Y-%m-%d %H:%M:%S"), "%Y-%m-%d"))

sampledata <- merge(sampledata, soils, all.x=T)

rownames(sampledata) <- sampledata$dnaSampleID

# Subset to samples present in both dataframes
common_samples <- intersect(rownames(sampledata), rownames(seqtab_collapse))

# Combine into phyloseq object
ps <- phyloseq(otu_table(seqtab_collapse[common_samples,], taxa_are_rows=FALSE),
               sample_data(sampledata[common_samples,]))

# store the DNA sequences of our ASVs in the refseq slot of the phyloseq object,
# and then rename our taxa to a short string
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Print phyloseq summary
ps


# Alpha diversity analysis
ps_rare <- rarefy_even_depth(ps, sample.size=10000, rngseed=101010)
alpha_div <- cbind("dnaSampleID" = get_variable(ps_rare, "dnaSampleID"), estimate_richness(ps_rare, measures = c("Observed", "Shannon")))
alpha_div <- merge(sampledata, alpha_div, by="dnaSampleID", all.y=TRUE)
alpha_div_long <- pivot_longer(alpha_div, c(Observed, Shannon), names_to="Metric", values_to="Value")
theme_set(theme_bw())
p1 <- gridExtra::arrangeGrob(
  ggplot(alpha_div, aes(x=siteID, y=Observed)) +
    geom_boxplot() +
    ylab("Observed ASV richness") +
    xlab("Site") +
    theme(panel.grid.major.x = element_blank()),
  ggplot(alpha_div, aes(x=siteID, y=Shannon)) +
    geom_boxplot() +
    ylab("ASV Shannon index") +
    xlab("Site") +
    theme(panel.grid.major.x = element_blank()),
  nrow=1
)
ggsave("./plots/grassland_alpha_diversity.png", plot=p1, device="png", width=4, height=2.5, units="in")

# Dissimilarity analysis
ps_prop <- transform_sample_counts(ps, function(x) x/sum(x))
set.seed(101010)
ps_ord <- ordinate(ps_prop, method = "NMDS")
p2 <- plot_ordination(ps_prop, ps_ord, color="soilInCaClpH", shape="siteID") +
  scale_color_gradient("Soil pH", low="red", high="blue") +
  scale_shape_manual("Site", values=c(21, 22, 24)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("./plots/grassland_ordination.png", plot=p2, device="png", width=5, height=4, units="in")
```
