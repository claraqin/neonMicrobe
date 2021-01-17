# plotEEvsTruncLen.R

library(ShortRead)
library(tidyr)
library(ggplot2)

source("R/utils.R") # load the plotEEProfile function
source("code/params.R")

runIDs <- c("runC5B2R", "runCBTWG", "runBFDG8")
# runIDs <-  c("runB69PP", "runB69RN", "runB9994", "runBDNB6", "runBF462", "runBFDG8", "runBNMJ5", "runBNMWB", "runBRPH4", "runC24VW", "runC25T6", "runC5B2R", "runC7WK3", "runC8VMV", "runC977L", "runC983L", "runCBJYB", "runCBTWG", "runCDHG2", "runCDJ5J")

if(is.null(PRESET_OUTDIR_SEQUENCE) | PRESET_OUTDIR_SEQUENCE == "") {
  PATH_16S <- file.path(PRESET_OUTDIR, "raw_sequence", "16S")
} else {
  PATH_16S <- file.path(PRESET_OUTDIR, PRESET_OUTDIR_SEQUENCE, "16S")
}
PATH_RAW <- file.path(PATH_16S, "0_raw")

PATH_OUTPUT <- file.path(PATH_16S, "qf_test_1-13-2021", "results")

rawFs <- rawRs <- list()
for(i in 1:length(runIDs)) {
  rawFs[[i]] <- list.files(PATH_RAW, pattern=paste0(runIDs[[i]], ".*_R1\\.fastq"), full.names=TRUE)
  rawRs[[i]] <- list.files(PATH_RAW, pattern=paste0(runIDs[[i]], ".*_R2\\.fastq"), full.names=TRUE)
}

# Use only first 10 samples from each run ID
rawFs <- lapply(rawFs, `[`, 1:10)
rawRs <- lapply(rawRs, `[`, 1:10)

g <- gridExtra::arrangeGrob(
  gridExtra::arrangeGrob(grobs = lapply(rawFs, function(x) plotEEProfile(x, aggregate=TRUE) + ggtitle(paste0(runIDs[i], ", R1")))),
  gridExtra::arrangeGrob(grobs = lapply(rawRs, function(x) plotEEProfile(x, aggregate=TRUE) + ggtitle(paste0(runIDs[i], ", R2")))),
  ncol=2
)

ggsave(file.path(PATH_OUTPUT, "expected_errors_plot.png"), plot=g, width=5, height=7, units="in")
