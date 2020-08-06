# Merge various NEON data products to link the outputs of the
# DADA2 pipeline (e.g. non-chimeric sequence table, taxonomic table)
# to environmental/sample data.
library(phyloseq)
library(ggplot2)
library(dplyr)

# Load parameters and utilities and create soil database if necessary
source("./code/utils.R")
source("./code/params.R")
#source("./code/zoey_params.R")

if(!file.exists(file.path(PRESET_OUTDIR_SOIL_DB, "soilDatabase.rds"))) { # Takes a while to run for the first time
  source("./code/create_soil_df.R")
} 
sampledata_full <- readRDS(file.path(PRESET_OUTDIR_SOIL_DB, "soilDatabase.rds"))

# Load combined sequence table and taxonomic table
seqtab.nochim <- readRDS("./data/NEON_16S_seqtab_nochim.Rds")
taxa <- readRDS("./data/NEON_16S_tax.Rds")

# Create table to link internalLabID with sampleID
link_geneticID_labID <- sampledata_full %>%
  select(geneticSampleID, internalLabID) %>%
  distinct() 

# Fix rownames on sequence table
rownames(seqtab.nochim) <- gsub("^run....._", "", rownames(seqtab.nochim))
rownames(seqtab.nochim) <- gsub("_16S_....._filt", "", rownames(seqtab.nochim))
sampledata <- sampledata_full[match(rownames(seqtab.nochim), sampledata_full$internalLabID),]
sampledata <- sampledata[, -which(names(sampledata)=="internalLabID")] 

# Remove rows without basic info
sampledata <- sampledata[which(!is.na(sampledata$sampleID)),]
seqtab.nochim <- seqtab.nochim[which(!is.na(sampledata$sampleID)),]

# Fix date column, add rownames for phyloseq
sampledata$date_ymd <- as.Date(format(as.Date(sampledata$collectDate, format="%Y-%m-%dT%H:%MZ"), "%Y-%m-%d"))
rownames(sampledata) <- sampledata$sampleID
rownames(seqtab.nochim) <- sampledata$sampleID

# Combine into phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(sampledata),
               tax_table(taxa))

# store the DNA sequences of our ASVs in the refseq slot of the phyloseq object,
# and then rename our taxa to a short string
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

saveRDS(ps, "./data/NEON_16S_phyloseq.rds")