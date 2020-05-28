# Create relational database using SQLite for the soil-related data
# in "./raw_data/sequence_metadata" and "./raw_data/soil"

library(RSQLite)

PRNUM_seqmeta <- "10108" # DP1.10108: "Soil microbe marker gene sequences"
PRNUM_chem <- "10078"    # DP1.10078: "Soil chemical properties (Distributed periodic)"
PRNUM_phys <- "10086"    # DP1.10086: "Soil physical properties (Distributed periodic)"
dpID_seqmeta <- paste("DP1", PRNUM_seqmeta, "001", sep=".")
dpID_chem <- paste("DP1", PRNUM_chem, "001", sep=".")
dpID_phys <- paste("DP1", PRNUM_phys, "001", sep=".")
stackDir_seqmeta <- file.path(PRESET_OUTDIR_SEQMETA, paste0("filesToStack",PRNUM_seqmeta))
stackDir_chem <- file.path(PRESET_OUTDIR_SOIL, paste0("filesToStack",PRNUM_chem))
stackDir_phys <- file.path(PRESET_OUTDIR_SOIL, paste0("filesToStack",PRNUM_phys))

downloadRawSoilData(return_data=FALSE, overwrite=TRUE)
# Warning: may take a while

path_rawDataFiles <- file.path(stackDir_seqmeta, "stackedFiles", "mmg_soilRawDataFiles.csv")
path_dnaExtraction <- file.path(stackDir_seqmeta, "stackedFiles", "mmg_soilDnaExtraction.csv")
path_soilCoreCollection <- paste(stackDir_phys, "stackedFiles", "sls_soilCoreCollection.csv", sep="/")
path_soilMoisture <- paste(stackDir_phys, "stackedFiles", "sls_soilMoisture.csv", sep="/")
path_soilpH <- paste(stackDir_phys, "stackedFiles", "sls_soilpH.csv", sep="/")
path_soilChemistry <- paste(stackDir_chem, "stackedFiles", "sls_soilChemistry.csv", sep="/")

dat_rawDataFiles <- read.delim(path_rawDataFiles, sep=",", stringsAsFactors=FALSE, na.strings = c("NA", ""))
dat_dnaExtraction <- read.delim(path_dnaExtraction, sep=",", stringsAsFactors=FALSE, na.strings = c("NA", ""))
dat_soilCoreCollection <- read.delim(path_soilCoreCollection, sep=",", stringsAsFactors=FALSE, na.strings = c("NA", ""))
dat_soilMoisture <- read.delim(path_soilMoisture, sep=",", stringsAsFactors=FALSE, na.strings = c("NA", ""))
dat_soilpH <- read.delim(path_soilpH, sep=",", stringsAsFactors=FALSE, na.strings = c("NA", ""))
dat_soilChemistry <- read.delim(path_soilChemistry, sep=",", stringsAsFactors=FALSE, na.strings = c("NA", ""))

# remove analytical replicates
dat_soilChemistry <- dat_soilChemistry[dat_soilChemistry$analyticalRepNumber == 1,]

# Fix unmerged %C and %N columns in soil moisture data
dup <- dat_soilChemistry[duplicated(dat_soilChemistry$sampleID),]
dup <- dat_soilChemistry[which(dat_soilChemistry$sampleID %in% dup$sampleID),] # all duplicated rows, including first appearances
dup %>% arrange(sampleID) %>%
  select(sampleID, analyticalRepNumber, nitrogenPercent, organicCPercent)
c.only.ind <- which(is.na(dup$nitrogenPercent) & !is.na(dup$organicCPercent))
n.only.ind <- which(!is.na(dup$nitrogenPercent) & is.na(dup$organicCPercent))
c_n_merged <- full_join(select(dup[c.only.ind,], -nitrogenPercent),
                        select(dup[n.only.ind,], sampleID, nitrogenPercent))

# merge back in with rest of C/N data
dat_soilChemistry <- dat_soilChemistry[!(dat_soilChemistry$sampleID %in% c_n_merged$sampleID),]
dat_soilChemistry$CNvals_merged <- FALSE
c_n_merged$CNvals_merged <- TRUE # so we know for later that these values are different
c_n_merged <- c_n_merged[,colnames(dat_soilChemistry)] # reorder column names
c_n_merged$CNratio <- round(c_n_merged$organicCPercent/c_n_merged$nitrogenPercent,1)
dat_soilChemistry <- rbind(dat_soilChemistry, c_n_merged) # merge back together

# Select a subset of cols from each table
joining_cols <- c("domainID", "siteID", "plotID", "sampleID")
dat_rawDataFiles <- select(
  dat_rawDataFiles, domainID, siteID, internalLabID, rawDataFileName, 
  rawDataFilePath, rawDataQF=dataQF
)
dat_dnaExtraction <- select(
  dat_dnaExtraction, domainID, siteID, plotID, geneticSampleID, internalLabID,
  extractionDataQF=dataQF
)
dat_soilCoreCollection <- select(
  dat_soilCoreCollection, domainID, siteID, plotID, namedLocation, plotType, nlcdClass, 
  coreCoordinateX, coreCoordinateY, geodeticDatum, decimalLatitude, decimalLongitude, 
  elevation, sccSamplingProtocolVersion=samplingProtocolVersion, collectDate, 
  sampleTiming, standingWaterDepth, nTransBoutType, sampleID, horizon, soilTemp, 
  litterDepth, sampleTopDepth, sampleBottomDepth, soilSamplingDevice, geneticSampleID, 
  sccDataQF=dataQF
)
dat_soilMoisture <- select(
  dat_soilMoisture, all_of(joining_cols), moistureSampleID, 
  smSamplingProtocolVersion=samplingProtocolVersion, soilMoisture, smDataQF
)
dat_soilpH <- select(
  dat_soilpH, all_of(joining_cols), pHSampleID, 
  pHSamplingProtocolVersion=samplingProtocolVersion, soilInWaterpH, soilInCaClpH, 
  pHDataQF
)
dat_soilChemistry <- select(
  dat_soilChemistry, all_of(joining_cols), cnSampleID, nitrogenPercent, 
  organicCPercent, CNratio, cnTestMethod=testMethod, cnInstrument=instrument, 
  cnDataQF=dataQF
)

# Add geolocation data to soilCoreCollection
library(geoNEON)
dat_soilCoreCollection <- getLocTOS(dat_soilCoreCollection, "sls_soilCoreCollection")

# Create database connection
con <- dbConnect(RSQLite::SQLite(), file.path(PRESET_OUTDIR_SOIL_DB, "soilDB.db"))
dbWriteTable(con, "rawDataFiles", dat_rawDataFiles, overwrite=TRUE)
dbWriteTable(con, "dnaExtraction", dat_dnaExtraction, overwrite=TRUE)
dbWriteTable(con, "soilCoreCollection", dat_soilCoreCollection, overwrite=TRUE)
dbWriteTable(con, "soilMoisture", dat_soilMoisture, overwrite=TRUE)
dbWriteTable(con, "soilpH", dat_soilpH, overwrite=TRUE)
dbWriteTable(con, "soilChemistry", dat_soilChemistry, overwrite=TRUE)
dbListTables(con)

