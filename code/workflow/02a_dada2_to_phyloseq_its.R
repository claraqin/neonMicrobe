# Merge various NEON data products to link the outputs of the
# DADA2 pipeline (e.g. non-chimeric sequence table, taxonomic table)
# to environmental/sample data.

# Load utilities and create soil database if necessary
source("./code/utils.R")
if(!file.exists(file.path(PRESET_OUTDIR_SOIL_DB, "soilDB.db"))) {
  source("./code/create_soil_database.R")
}

# Load libraries
library(phyloseq)
library(ggplot2)
library(RSQLite)

# Load combined sequence table and taxonomic table
seqtab.nochim <- readRDS(file.path(PRESET_OUTDIR_DADA2, PRESET_FILENAME_JOINED_SEQTAB)) # TODO: Modify these to be the files that you're interested in processing.
taxa <- readRDS(file.path(PRESET_OUTDIR_DADA2, PRESET_FILENAME_TAXTAB))

# Connect to soil database
con <- dbConnect(RSQLite::SQLite(), file.path(PRESET_OUTDIR_SOIL_DB, "soilDB.db"))
res <- dbSendQuery(con, "
  SELECT internalLabID, scc.namedLocation, scc.domainID, scc.siteID, scc.plotID,
         scc. plotType, scc.nlcdClass, scc.collectDate, scc.sampleTiming,
         scc.standingWaterDepth, scc.soilTemp, scc.horizon,
         scc.sampleID, sph.soilInCaClpH, sm.soilMoisture, sch.nitrogenPercent,
         sch.organicCPercent, sch.CNratio,
         scc.utmZone, scc.adjNorthing, scc.adjEasting, scc.adjDecimalLatitude,
         scc.adjDecimalLongitude, scc.adjElevation
  FROM dnaExtraction AS dna
  INNER JOIN soilCoreCollection AS scc
     ON dna.geneticSampleID=scc.geneticSampleID
  LEFT JOIN soilpH AS sph
     ON sph.sampleID=scc.sampleID
  LEFT JOIN soilMoisture AS sm
     ON sm.sampleID=scc.sampleID
  LEFT JOIN soilChemistry AS sch
     ON sch.sampleID=scc.sampleID
")
res_df <- dbFetch(res)
apply(res_df, 2, function(x){mean(is.na(x))})
dbClearResult(res)
dbDisconnect(con)

sampledata <- res_df[match(rownames(seqtab.nochim), res_df$internalLabID),]
sampledata <- sampledata[, -which(names(sampledata)=="internalLabID")]
rownames(sampledata) <- rownames(seqtab.nochim)

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

saveRDS(ps, "./data/NEON_ITS_phyloseq_DL08-13-2019.Rds")

# # Phyloseq object is ready for analysis
# 
# plot_richness(ps, x="collectDate", measures=c("Shannon", "Simpson"), color="siteID")
# 
# plot_richness(ps, x="soilTemp", measures=c("Shannon", "Simpson"), color="horizon")
# 
# plot_richness(ps, x="decimalLatitude", measures=c("Observed", "Shannon"), color="domainID") + 
#   geom_smooth()
# 
# # Subset for 2018 data
# dates <- get_variable(ps, "collectDate.x")
# dates_yr <- as.integer(format(as.Date(dates, format="%Y-%m-%dT%H:%MZ"), "%Y"))
# dates_mo <- as.integer(format(as.Date(dates, format="%Y-%m-%dT%H:%MZ"), "%m"))
# sample_data(ps)$year <- dates_yr
# sample_data(ps)$mo <- dates_mo
# ps_2017 <- subset_samples(ps, year == "2017")
# rm(ps_2017)
# ps_2017_06 <- subset_samples(ps, year == "2017" & mo == "6")
# 
# saveRDS(ps_2017_06, "./data/NEON_ITS_phyloseq_subset_2017_06.Rds")
