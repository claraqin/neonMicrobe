# Code to generate ITS sequence metadata used in examples

library(neonMicrobe)
library(neonUtilities)
library(dplyr)
library(usethis)

seqmeta_greatplains_its <- downloadSequenceMetadata(startYrMo = "2017-07", endYrMo = "2017-07",
                                     sites = c("KONZ", "CPER", "NOGP"),
                                     targetGene = "ITS", outDir=FALSE)
write.csv(seqmeta_greatplains_its, "data-raw/seqmeta_greatplains_its.csv")
write.csv(seqmeta_greatplains_its, "inst/extdata/seqmeta_greatplains_its.csv")
usethis::use_data(seqmeta_greatplains_its, overwrite=TRUE)
