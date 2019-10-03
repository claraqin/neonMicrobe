# Merge various NEON data products to link the outputs of the
# DADA2 pipeline (e.g. non-chimeric sequence table, taxonomic table)
# to environmental/sample data.

source("./utils.R")
library(phyloseq)
library(ggplot2)

seqtab.nochim <- readRDS("./NEON_ITS_seqtab_nochim_DL08-13-2019.Rds")
taxa <- readRDS("./NEON_ITS_taxa_DL08-13-2019.Rds")

seqmetadata <- downloadAllSequenceMetadata()

seqmetadata$geneticSampleID <- sub("-DNA[1-3]", "", seqmetadata$dnaSampleID)

# seqmetadata %>%
#   select(geneticSampleID, internalLabID) %>%
#   distinct() ->
#   link_geneticID_labID

seqmetadata %>%
  select(geneticSampleID, rawDataFileName) %>%
  mutate(filename = sub(".fastq"))
  distinct() ->
  link_geneticID_filename

# TODO: Merge metadata with the output of downloadAllRawSoilData(),
# joining by geneticSampleID. Then match rownames of seqtab_nochim
# to the environmental data via the newly merged dataset's 
# "internalLabID" column.

soildata <- downloadAllRawSoilData()

# soildata %>%
#   filter(geneticSampleID != "") %>%
#   inner_join(seqmetadata, by="geneticSampleID") ->
#   sampledata

soildata %>%
  filter(geneticSampleID != "") %>%
  full_join(link_geneticID_labID, by="geneticSampleID") ->
  sampledata_full

seqtab_labID <- rownames(seqtab.nochim)

sampledata_ind <- match(seqtab_labID, sampledata_full$internalLabID)
# Not all sequencing samples have a corresponding sample-data entry!

sampledata <- sampledata_full[sampledata_ind,]
rownames(sampledata) <- seqtab_labID


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

saveRDS(ps, "NEON_ITS_phyloseq_DL08-13-2019.Rds")

# Phyloseq object is ready for analysis

plot_richness(ps, x="collectDate", measures=c("Shannon", "Simpson"), color="siteID")

plot_richness(ps, x="soilTemp", measures=c("Shannon", "Simpson"), color="horizon")

plot_richness(ps, x="decimalLatitude", measures=c("Observed", "Shannon"), color="domainID")

