# plotEEvsTruncLen.R

library(ShortRead)
library(tidyr)
library(ggplot2)

source("R/utils.R") # load the plotEEProfile function
source("code/params.R")

runIDs <- c("runC5B2R", "runCBTWG", "runBFDG8")

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

g <- gridExtra::arrangeGrob(
  gridExtra::arrangeGrob(
    plotEEProfile(rawFs[[1]][1:10], aggregate=TRUE) + ggtitle("C5B2R, R1"),
    plotEEProfile(rawFs[[2]][1:10], aggregate=TRUE) + ggtitle("CBTWG, R1"),
    plotEEProfile(rawFs[[3]][1:10], aggregate=TRUE) + ggtitle("BFDG8, R1")
  ),
  gridExtra::arrangeGrob(
    plotEEProfile(rawRs[[1]][1:10], aggregate=TRUE) + ggtitle("C5B2R, R2"),
    plotEEProfile(rawRs[[2]][1:10], aggregate=TRUE) + ggtitle("CBTWG, R2"),
    plotEEProfile(rawRs[[3]][1:10], aggregate=TRUE) + ggtitle("BFDG8, R2")
  ),
  ncol=2
)

ggsave(file.path(PATH_OUTPUT, "expected_errors_plot.png"), plot=g, width=5, height=7, units="in")
