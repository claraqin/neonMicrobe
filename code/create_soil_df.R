# Create relational database using SQLite for the soil-related data
# in "./raw_data/sequence_metadata" and "./raw_data/soil"
library(geoNEON)

# Download NEON's soil data and sequence metadata
# Takes ~5 minutes to run if not previously downloaded
dat_soil <- downloadRawSoilData(sites = 'all', startYrMo = "2013-06", endYrMo = "2020-07", return_data=TRUE, overwrite=FALSE)

# Takes ~4 minutes to run if not previously downloaded
# TO-DO: The output of downloadSequenceMetadataRev() includes download date - which makes it harder to read in flexibly 
path_seqMeta <- file.path(PRESET_OUTDIR_SEQMETA, "mmg_soilMetadata_ITS_2020-08-06.csv")
if (file.exists(path_seqMeta)){
  dat_sequenceMetadata <- read.csv(path_seqMeta, stringsAsFactors=FALSE, na.strings = c("NA", ""))
  } else {
  dat_sequenceMetadata <- downloadSequenceMetadataRev(sites = 'all', startYrMo = "2013-06", endYrMo = "2020-07", 
                                                    targetGene = "ITS", dir = PRESET_OUTDIR_SEQMETA)
}

# remove analytical replicates
dat_soil <- dat_soil[which(!dat_soil$analyticalRepNumber %in% c(2:4)),]

# Fix unmerged %C and %N columns in soil moisture data
dup <- dat_soil[duplicated(dat_soil$cnSampleID) & !is.na(dat_soil$cnSampleID),]
dup <- dat_soil[which(dat_soil$cnSampleID %in% dup$cnSampleID),] # all duplicated rows, including first appearances
c.only.ind <- which(is.na(dup$nitrogenPercent) & !is.na(dup$organicCPercent))
n.only.ind <- which(!is.na(dup$nitrogenPercent) & is.na(dup$organicCPercent))
c_n_merged <- full_join(select(dup[c.only.ind,], -nitrogenPercent),
                        select(dup[n.only.ind,], sampleID, nitrogenPercent))
c_n_merged <- c_n_merged[,colnames(dat_soil)] # reorder column names
c_n_merged$CNratio <- round(c_n_merged$organicCPercent/c_n_merged$nitrogenPercent,1)
c_n_merged$CNvals_merged <- TRUE # so we know for later that these values are different

# merge back in with rest of C/N data
dat_soil$CNvals_merged <- FALSE
dat_soil <- dat_soil[!(dat_soil$sampleID %in% c_n_merged$sampleID),]
dat_soil <- rbind(dat_soil, c_n_merged)

# Add geolocation data to soilCoreCollection
dat_soil <- getLocTOS(dat_soil, "sls_soilCoreCollection")

# Select a subset of cols from sequence metadata
dat_sequenceMetadata <- select(
  dat_sequenceMetadata, domainID = domainID.x, siteID = siteID.x, plotID, dnaSampleID, 
  rawDataFileName, rawDataFilePath, internalLabID,
  extractionDataQF=dataQF.x, rawDataQF=dataQF.y
)
dat_sequenceMetadata$geneticSampleID <- gsub("-DNA.$", "", dat_sequenceMetadata$dnaSampleID)

# Merge: only keeping soil data if there is accompanying sequence data (won't work for legacy sequence data)
common_cols <- intersect(colnames(dat_sequenceMetadata), colnames(dat_soil))
out <- merge(dat_sequenceMetadata, dat_soil, all.x = T, by = common_cols) 
#out2 <- merge(dat_sequenceMetadata, dat_soil, all = T, by = common_cols)

# Remove duplicates and select output columns 
out <- out[!duplicated(out),]
out <- out %>% select(domainID, siteID, plotID, namedLocation, plotType, nlcdClass, adjDecimalLatitude, adjDecimalLongitude, 
  adjElevation, collectDate, sampleTiming, standingWaterDepth, sampleID, horizon, soilTemp, 
  litterDepth, sampleTopDepth, sampleBottomDepth, geneticSampleID,nitrogenPercent, organicCPercent, CNratio, CNvals_merged, soilInWaterpH, soilInCaClpH
)

saveRDS(out, file = file.path(PRESET_OUTDIR_SOIL_DB, "soilDatabase.rds"))

