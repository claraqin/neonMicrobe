# Preset parameters for utils.R, new_server_setup.R, dada2_workflow_its.R, dada2_workflow_16s.R, and reorg_sequence_data.R

## PARAMETERS FOR NEW_SERVER_SETUP
PRESET_SITES = "all"
PRESET_START_YR_MO = "2011-06"
PRESET_END_YR_MO = format(Sys.Date(), "%Y-%m")
PRESET_OUTDIR_SEQUENCE = "/data/ZHULAB/NEON_DOB/Illumina/NEON"
PRESET_OUTDIR_SEQMETA = "/data/ZHULAB/NEON_DOB/sequence_metadata"
PRESET_OUTDIR_SOIL = "/data/ZHULAB/NEON_DOB/soil"
PRESET_CHECK_FILE_SIZE = FALSE
PRESET_RETURN_DATA = TRUE
SITE_AND_DATE_RANGE_FILENAME = "sites_and_date_range.txt"

## PARAMETERS FOR DADA2_WORKFLOW
# Primers
PRIMER_ITS_FWD = "CTTGGTCATTTAGAGGAAGTAA" # Forward primer sequence
PRIMER_ITS_REV = "GCTGCGTTCTTCATCGATGC" # Reverse primer sequence
# Cutadapt path
CUTADAPT_PATH = "/afs/cats.ucsc.edu/users/b/claraqin/.local/bin/cutadapt"
# Whether to download and process only a small subset
SMALL_SUBSET = FALSE # If TRUE, all fastq.gz.tar files matching these parameters will still be downloaded,
                     # and all fastq files within them will be unzipped, but only the first 2 forward-reverse
                     # pairs of sequence files from the first sequencing run ID will be processed in 
                     # dada2_workflow.R and in following scripts.
VERBOSE = FALSE # If TRUE, will generate additional output in dada2_workflow_its.R
# DADA2
