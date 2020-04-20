# Merge various NEON data products to link the outputs of the
# DADA2 pipeline (e.g. non-chimeric sequence table, taxonomic table)
# to environmental/sample data.

source("./code/utils.R")
library(phyloseq)
library(ggplot2)

seqtab.nochim <- readRDS(file.path(PRESET_OUTDIR_DADA2, PRESET_FILENAME_JOINED_SEQTAB)) # TODO: Modify these to be the files that you're interested in processing.
taxa <- readRDS(file.path(PRESET_OUTDIR_DADA2, PRESET_FILENAME_TAXTAB))

# seqmetadata <- downloadSequenceMetadata()
# 
# seqmetadata$geneticSampleID <- sub("-DNA[1-3]", "", seqmetadata$dnaSampleID)
# 
# seqmetadata %>%
#   # select(geneticSampleID, rawDataFileName, internalLabID) %>%
#   # mutate(match_name = sub("_((ITS)|(16S))_R[12].fastq(.tar)?(.gz)?$","", rawDataFileName)) %>% dim()
#   select(geneticSampleID, internalLabID) %>%
#   distinct() ->
#   link_geneticID_labID

soildata <- downloadRawSoilData()

# soildata %>%
#   filter(geneticSampleID != "") %>%
#   full_join(link_geneticID_labID, by="geneticSampleID") ->
#   sampledata_full

seqtab_internalLabID <- rownames(seqtab.nochim)

# sampledata_ind <- match(seqtab_internalLabID, sampledata_full$internalLabID)
# # Not all sequencing samples have a corresponding sample-data entry!

# Try again with Lee's custom metadata (TODO: need to operationalize)
lee_metadata <- read.csv("./code/ITSmetadataMapped.csv")
lee_metadata %>%
  mutate(match_name = sub("_((ITS)|(16S))_R[12].fastq$", "", basename(as.character(fileName)))) ->
  lee_metadata
sampledata_ind2 <- match(seqtab_internalLabID, lee_metadata$match_name)
mean(is.na(sampledata_ind2)) # 0.262
# Yes, let's go with this one instead of seqmetadata

soildata %>%
  filter(geneticSampleID != "") %>%
  full_join(select(lee_metadata, match_name, geneticSampleID), by="geneticSampleID") ->
  sampledata_full
sampledata_ind <- match(seqtab_internalLabID, sampledata_full$match_name)

sampledata <- sampledata_full[sampledata_ind,]
rownames(sampledata) <- seqtab_internalLabID


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

saveRDS(ps, "./data/NEON_ITS_phyloseq_DL08-13-2019.Rds")

# Phyloseq object is ready for analysis

plot_richness(ps, x="collectDate", measures=c("Shannon", "Simpson"), color="siteID")

plot_richness(ps, x="soilTemp", measures=c("Shannon", "Simpson"), color="horizon")

plot_richness(ps, x="decimalLatitude", measures=c("Observed", "Shannon"), color="domainID") + 
  geom_smooth()

# Subset for 2018 data
dates <- get_variable(ps, "collectDate.x")
dates_yr <- as.integer(format(as.Date(dates, format="%Y-%m-%dT%H:%MZ"), "%Y"))
dates_mo <- as.integer(format(as.Date(dates, format="%Y-%m-%dT%H:%MZ"), "%m"))
sample_data(ps)$year <- dates_yr
sample_data(ps)$mo <- dates_mo
ps_2017 <- subset_samples(ps, year == "2017")
rm(ps_2017)
ps_2017_06 <- subset_samples(ps, year == "2017" & mo == "6")

saveRDS(ps_2017_06, "./data/NEON_ITS_phyloseq_subset_2017_06.Rds")




write.csv(table(ps_sampledata$siteID, ps_sampledata$year), "./data/sites_by_year.csv")

table(sample_data(ps_2017_06)$siteID)

write.csv(table(ps_sampledata$year, ps_sampledata$mo), "./data/obs_by_year_and_mo.csv")
