# Merge various NEON data products to link the outputs of the
# DADA2 pipeline (e.g. non-chimeric sequence table, taxonomic table)
# to environmental/sample data.

source("./code/utils.R")
library(phyloseq)
library(ggplot2)
library(dplyr)

seqtab.nochim <- readRDS("./data/NEON_16S_seqtab_nochim_DL08-13-2019.Rds")
taxa <- readRDS("./data/NEON_16S_taxa_DL08-13-2019.Rds")

seqmetadata <- downloadAllSequenceMetadata()

seqmetadata$geneticSampleID <- sub("-DNA[1-3]", "", seqmetadata$dnaSampleID)

seqmetadata %>%
  # select(geneticSampleID, rawDataFileName, internalLabID) %>%
  # mutate(match_name = sub("_((ITS)|(16S))_R[12].fastq(.tar)?(.gz)?$","", rawDataFileName)) %>% dim()
  select(geneticSampleID, internalLabID) %>%
  distinct() ->
  link_geneticID_labID

soildata <- downloadAllRawSoilData()

soildata %>%
  filter(geneticSampleID != "") %>%
  full_join(link_geneticID_labID, by="geneticSampleID") ->
  sampledata_full

seqtab_rownames <- sub("_16S_(F|R)_filt.fastq.gz$", "", rownames(seqtab.nochim))

sampledata_ind <- match(seqtab_rownames, sampledata_full$internalLabID)
# # Not all sequencing samples have a corresponding sample-data entry!

# # Try again with Lee's custom metadata (TODO: need to operationalize)
# lee_metadata <- read.csv("./code/ITSmetadataMapped.csv")
# lee_metadata %>%
#   mutate(match_name = sub("_((ITS)|(16S))_R[12].fastq$", "", basename(as.character(fileName)))) ->
#   lee_metadata
# sampledata_ind2 <- match(seqtab_rownames, lee_metadata$match_name)
# mean(is.na(sampledata_ind2)) # 0.262
# # Yes, let's go with this one instead of seqmetadata
#
# soildata %>%
#   filter(geneticSampleID != "") -> # %>%
#   # full_join(select(lee_metadata, match_name, geneticSampleID), by="geneticSampleID") ->
#   sampledata_full
# sampledata_ind <- match(seqtab_rownames, sampledata_full$match_name)

sampledata <- sampledata_full[sampledata_ind,]

sampledata %>%
  select(sampleID,siteID,plotID, collectDate.x) %>%
  # mutate(collectDate.x = as.Date(collectDate.x, "%Y-%m-%dT%H:%MZ")) %>%
  mutate(sampleID = as.factor(sampleID),
         siteID = as.factor(siteID),
         plotID = as.factor(plotID)) ->
  sampledata_select
rownames(sampledata_select) <- rownames(seqtab.nochim)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(sampledata_select),
               tax_table(taxa))

# store the DNA sequences of our ASVs in the refseq slot of the phyloseq object,
# and then rename our taxa to a short string
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

saveRDS(ps, "./data/NEON_16S_phyloseq_DL08-13-2019.Rds")

# ps <- readRDS("./data/NEON_16S_phyloseq_DL08-13-2019.Rds")

ps_dated <- subset_samples(ps, !is.na(collectDate.x))
dates <- get_variable(ps_dated, "collectDate.x")
dates_ymd <- format(as.Date(dates, format="%Y-%m-%dT%H:%MZ"), "%Y-%m-%d")
sample_data(ps_dated)$date_ymd <- dates_ymd

ps_genus <- tax_glom(ps_dated, taxrank = "Genus")

ps_genus_rm0 <- prune_samples(sample_sums(ps_genus) > 0, ps_genus)

# ps_genus_aggsample <- merge_samples(ps_genus, "plotID")

# ps_genus_prop <- transform_sample_counts(ps_genus_aggsample, function(OTU) OTU/sum(OTU))
ps_genus_prop <- transform_sample_counts(ps_genus_rm0, function(OTU) OTU/sum(OTU))

y <- as(otu_table(ps_genus_prop), "matrix")
hist(y)
hist(y[y!=0])
y_nonzero_per_sample <- apply(y, 1, function(x) { sum(x != 0, na.rm=TRUE) })
hist(y_nonzero_per_sample)
mean(y_nonzero_per_sample > 1) # 65% of all samples contain more than one genus

ps_genus_melt <- psmelt(ps_genus_prop)
# ps_genus_melt <- psmelt(ps_genus)

ps_genus_melt %>%
  select(sampleID,siteID,plotID,date_ymd,Genus,Abundance) ->
  out
head(out)
mean(is.na(out$Abundance))
mean(out$Abundance == 0)
hist(out$Abundance[out$Abundance != 0])

out_rmna <- out[!is.na(out$Abundance),]
write.csv(out_rmna, file = "./data/NEON_16S_dated_genus_abundance.csv", row.names=FALSE)

asv_richness <- apply(otu_table(ps_dated),1, function(x){ sum(x!=0) })
hist(asv_richness)
head(asv_richness)
asv_richness_df <- data.frame(
  sampleID = get_variable(ps_dated, "sampleID"),
  siteID = get_variable(ps_dated, "siteID"),
  plotID = get_variable(ps_dated, "plotID"),
  date_ymd = get_variable(ps_dated, "date_ymd"),
  asv_richness = asv_richness
)
head(asv_richness_df)
write.csv(asv_richness_df, file = "./data/NEON_16S_asv_richness.csv", row.names=FALSE)
