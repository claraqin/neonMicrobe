# Code to generate sequence metadata used in examples

library(neonMicrobe)
library(neonUtilities)
library(dplyr)
library(usethis)

seqmeta_greatplains <- downloadSequenceMetadata(startYrMo = "2017-07", endYrMo = "2017-07",
                                     sites = c("KONZ", "CPER", "NOGP"),
                                     targetGene = "16S", outDir=FALSE)
write.csv(seqmeta_greatplains, "data-raw/seqmeta_greatplains.csv")
write.csv(seqmeta_greatplains, "inst/extdata/seqmeta_greatplains.csv")
usethis::use_data(seqmeta_greatplains, overwrite=TRUE)
