# Preset parameters for utils.R, new_server_setup.R, dada2_workflow_its.R, dada2_workflow_16s.R, and reorg_sequence_data.R

## PARAMETERS FOR NEW_SERVER_SETUP AND OUTPUT DIRECTORIES

# Output directory for the pipeline. Absolute path.
# Output directory for the pipeline. Creates a "NEON" directory within current working directory.
PRESET_OUTDIR = file.path(getwd(), "NEON")

# Can leave the following blank to generate default directory structure.
# Or, character string to append to the end of PRESET_OUTDIR to generate
# and use alternative directory structures.
PRESET_OUTDIR_SEQUENCE = "raw_sequence" # for sequence data (fastq files)
PRESET_OUTDIR_SEQMETA = "sequence_metadata" # for sequence metadata
PRESET_OUTDIR_SOIL = "soil" # for soil data

# Other output directories and filenames. Absolute paths.
PRESET_OUTDIR_SOIL_DB = "/data/ZHULAB/NEON_DOB" # Database containing both seqmeta and soil data
PRESET_OUTDIR_DADA2 = "/raid/users/claraqin/zhulab/neonSoilMicrobeProcessing/data" # for phyloseq outputs
PRESET_FILENAME_JOINED_SEQTAB = "NEON_ITS_seqtab_nochim_DL08-13-2019_truncQ4_maxEE8_R1only.Rds"
PRESET_FILENAME_TAXTAB = "NEON_ITS_taxa_DL08-13-2019_truncQ4_maxEE8_R1only.Rds"

## PARAMETERS FOR DADA2_WORKFLOW

# Primers
PRIMER_ITS_FWD = "CTTGGTCATTTAGAGGAAGTAA" # Forward primer sequence
PRIMER_ITS_REV = "GCTGCGTTCTTCATCGATGC" # Reverse primer sequence
PRIMER_16S_FWD = "CCTACGGGNBGCASCAG" # Forward primer sequence
PRIMER_16S_REV = "GACTACNVGGGTATCTAATCC" # Reverse primer sequence

# Search for & set the Cutadapt path
# system2 should work on Windows, but I'm not sure if the "which" command works.
CUTADAPT_PATH <- system2("which", args = "cutadapt", stdout = TRUE)
if(length(CUTADAPT_PATH) == 0){
  message("Could not find the cutadapt tool in file system. This tool is necessary for ITS (fungal) raw sequence processing. Please load the module if necessary (i.e. 'module load cutadapt') or download from https://cutadapt.readthedocs.io/en/stable/installation.html")
}

# UNITE reference database (FASTA file) path
UNITE_REF_PATH = "/data/sh_general_release_dynamic_04.02.2020.fasta"
SILVA_REF_PATH = "/projectnb/talbot-lab-data/zrwerbin/decomposition/silva_training_set.fa.gz"

# Whether to download and process only a small subset
SMALL_SUBSET = FALSE # If TRUE, all fastq.gz.tar files matching these parameters will still be downloaded,
# and all fastq files within them will be unzipped, but only the first 2 forward-reverse
# pairs of sequence files from the first sequencing run ID will be processed in
# dada2_workflow.R and in following scripts.

# Whether to generate additional output in dada2_workflow_its.R
VERBOSE = FALSE

# Whether to use multithreading in DADA2 (Windows users should set to FALSE),
# or, if integer is provided, how many threads to use
#MULTITHREAD = TRUE
if (Sys.info()["sysname"] == "Windows"){
  MULTITHREAD = FALSE
} else {
  MULTITHREAD = 4
}

# filterAndTrim arguments
MAX_EE_FWD = 8 # max. allowable expected errors in forward reads that pass filter
MAX_EE_REV = 8 # max. allowable expected errors in reverse reads that pass filter
TRUNC_Q = 2 # base quality score after which to truncate sequence
# NOTE: it may be desirable to set this higher than 2 if the goal
# is to increase the proportion of reads passing the filter, as this
# would allow fewer "expected errors" to appear in the sequence
# (by virtue of having shorter reads)

MIN_LEN = 50 # min. allowable length of reads that pass filter
# NOTE: it may be desirable to set this higher to increase the
# likelihood of sufficient overlap between read pairs. However,
# this is at odds with the incentive for setting TRUNC_Q higher
# (see previous)

# VALIDITY CHECKS (do not modify)
# if(!(TARGET_GENE %in% c("ITS", "16S", "all"))) warning("TARGET_GENE must be 'ITS', '16S', or 'all'")
