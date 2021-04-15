# process-its-sequences-cl.R

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  stop("No arguments provided. Must enter args as 'Rscript process-its-sequences-cl.R [arg1]', ",
       "where [arg1] is index of NEON sequencing run ID to be using, 1 to 20.")
}

library(plyr)
library(dplyr)
library(tibble)
library(tidyr)
library(Biostrings)
library(ShortRead)
library(phyloseq)
library(dada2)
library(ggplot2)
library(neonMicrobe)

print(getwd())
setBaseDirectory(getwd())
makeDataDirectories(check_location=FALSE)

i <- as.numeric(args[1])

meta <- read.csv("data/sequence_metadata/qc_metadata/mmg_metadata_ITS_QCd_20210331.csv",
                 stringsAsFactors = FALSE)

unique_runs <- sort(unique(meta$sequencerRunID))
runID <- unique_runs[i]
message("Processing run: ", runID)

meta_onerun <- meta[which(meta$sequencerRunID==runID),]
nrow(meta_onerun)

fl_nm <- meta_onerun$rawDataFileName
head(fl_nm)

if(runID == "BNM6G") {
  fl_nm <- fl_nm[-which(fl_nm == "BMI_Plate37WellA10_ITS_BNM6G_R1.fastq.gz")]
}

prefilter_trackReads <- qualityFilterITS(
  fl_nm, "raw", "1_filtN", meta = meta_onerun, maxN = 0, 
  multithread = TRUE # set multithread = FALSE on Windows computers though
)

trimPrimersITS(
  fl_nm, in_subdir = "1_filtN", out_subdir = "2_trimmed_once", 
  meta = meta_onerun, cutadapt_path = "/hb/home/claraqin/.conda/envs/cutadaptenv/bin/cutadapt"
)
trimPrimersITS(
  fl_nm, in_subdir = "2_trimmed_once", out_subdir = "2_trimmed", 
  meta = meta_onerun, cutadapt_path = "/hb/home/claraqin/.conda/envs/cutadaptenv/bin/cutadapt",
  primer_ITS_fwd = "GCTGCGTTCTTCATCGATGC", # swapping the primers
  primer_ITS_rev = "CTTGGTCATTTAGAGGAAGTAA"
)

filter_trackReads <- qualityFilterITS(
  fl_nm, in_subdir = "2_trimmed", out_subdir = "3_filtered",
  meta = meta_onerun, multithread = TRUE,
  maxN = 0, truncQ = 2, minLen = 50, maxEE = 8
)

dada_out <- runDadaITS(
  fl_nm, in_subdir = "3_filtered", meta = meta_onerun, 
  multithread = TRUE, remove_chimeras = TRUE,
  out_seqtab = file.path(NEONMICROBE_DIR_MIDPROCESS(), "ITS", "4_seqtabs", 
                         paste0("asv_its_", runID, ".Rds")),
  out_track = FALSE
)

track <- combineReadTrackingTablesITS(prefilter_trackReads, filter_trackReads, dada_out$track)
write.csv(track, file.path(NEONMICROBE_DIR_TRACKREADS(), "ITS",
                           paste0("track_its_", runID, ".csv")))

