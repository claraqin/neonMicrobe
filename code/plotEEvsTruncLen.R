# plotEEvsTruncLen.R

library(ShortRead)
library(tidyr)
library(ggplot2)
library(dplyr)

devtools::load_all() # TODO: Replace with library(neonMicrobe) after publishing

runIDs <-  c("B69PP", "B69RN", "B9994", "BDNB6", "BF462", "BFDG8", "BNMJ5", "BNMWB", "BRPH4", "C24VW", "C25T6", "C5B2R", "C7WK3", "C8VMV", "C977L", "C983L", "CBJYB", "CBTWG", "CDHG2", "CDJ5J")

rawFs <- rawRs <- list()
for(i in 1:length(runIDs)) {
  rawFs[[i]] <- list.files(file.path(NEONMICROBE_DIR_SEQUENCE(), "16S"), pattern=paste0(runIDs[[i]], ".*_R1\\.fastq"), full.names=TRUE)
  rawRs[[i]] <- list.files(file.path(NEONMICROBE_DIR_SEQUENCE(), "16S"), pattern=paste0(runIDs[[i]], ".*_R2\\.fastq"), full.names=TRUE)
}

# Use only first 10 samples from each run ID
rawFs <- lapply(rawFs, `[`, 1:10)
rawRs <- lapply(rawRs, `[`, 1:10)

theme_set(theme_bw())
g <- gridExtra::arrangeGrob(
  gridExtra::arrangeGrob(grobs = lapply(seq_along(runIDs), function(i) plotEEProfile(rawFs[[i]], aggregate=TRUE) + ggtitle(paste0(runIDs[i], ", R1"))), ncol=1),
  gridExtra::arrangeGrob(grobs = lapply(seq_along(runIDs), function(i) plotEEProfile(rawRs[[i]], aggregate=TRUE) + ggtitle(paste0(runIDs[i], ", R2"))), ncol=1),
  ncol=2
)

saveRDS(g, file.path(PATH_OUTPUT, "eeplot.Rds"))
# ggsave(file.path(PATH_OUTPUT, "expected_errors_plot.png"), plot=g, width=5, height=42, units="in")

