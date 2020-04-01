# Preset parameters for utils.R, new_server_setup.R, dada2_workflow_its.R, dada2_workflow_16s.R, and reorg_sequence_data.R

## PARAMETERS FOR NEW_SERVER_SETUP AND OUTPUTP DIRECTORIES
PRESET_SITES = "all"
PRESET_START_YR_MO = "2011-06"
PRESET_END_YR_MO = format(Sys.Date(), "%Y-%m")
PRESET_OUTDIR_SEQUENCE = "/data/ZHULAB/NEON_DOB/Illumina/NEON"
PRESET_OUTDIR_SEQMETA = "/data/ZHULAB/NEON_DOB/sequence_metadata"
PRESET_OUTDIR_SOIL = "/data/ZHULAB/NEON_DOB/soil"
PRESET_OUTDIR_DADA2 = "/raid/users/claraqin/zhulab/NEON_soil_microbe_processing/data"
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
# Whether to use multithreading in DADA2 (Windows users should set to FALSE),
# or, if integer is provided, how many threads to use
MULTITHREAD = FALSE
# filterAndTrim arguments
MAX_EE_FWD = 2 # max. allowable expected errors in forward reads that pass filter
MAX_EE_REV = 2 # max. allowable expected errors in reverse reads that pass filter
TRUNC_Q = 2 # base quality score after which to truncate sequence
            # NOTE: it may be desirable to set this higher than 2 if the goal
            # is to increase the proportion of reads passing the filter, as this
            # would allow fewer "expected errors" to appear in the sequence
            # (by virtue of having shorter sequences)
MIN_LEN = 50 # min. allowable length of sequences that pass filter
             # NOTE: it may be desirable to set this higher to increase the
             # likelihood of sufficient overlap between read pairs. However,
             # this is at odds with the incentive for setting TRUNC_Q higher
             # (see previous)
