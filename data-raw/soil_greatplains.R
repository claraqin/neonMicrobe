library(plyr)
library(dplyr)
library(neonUtilities)
library(R.utils)
library(neonMicrobe)
library(tidyr)

soil_greatplains <- downloadSoilData(startYrMo = "2017-07", endYrMo = "2017-07",
                                     sites = c("KONZ", "CPER", "NOGP"))
write.csv(seqmeta_greatplains, "data-raw/soil_greatplains.csv")
usethis::use_data(soil_greatplains, overwrite=TRUE)
