# Code to generate 16S sequence metadata used in examples

library(neonMicrobe)
library(neonUtilities)
library(dplyr)
library(usethis)

seqmeta_greatplains_16s <- downloadSequenceMetadata(startYrMo = "2017-07", endYrMo = "2017-07",
                                     sites = c("KONZ", "CPER", "NOGP"),
                                     targetGene = "16S", outDir=FALSE)
write.csv(seqmeta_greatplains_16s, "data-raw/seqmeta_greatplains_16s.csv")
write.csv(seqmeta_greatplains_16s, "inst/extdata/seqmeta_greatplains_16s.csv")
usethis::use_data(seqmeta_greatplains_16s, overwrite=TRUE)
