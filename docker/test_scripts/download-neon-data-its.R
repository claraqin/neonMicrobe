
library(plyr)
library(dplyr)
library(neonUtilities)
library(neonMicrobe)

print(getwd())
setBaseDirectory(getwd())
makeDataDirectories(check_location=FALSE)

meta_its <- downloadSequenceMetadata(targetGene = "ITS")

meta_its_qc <- qcMetadata(meta_its, pairedReads = "N", rmFlagged = "Y")

downloadRawSequenceData(meta_its_qc, verbose = TRUE, checkSize=FALSE)
