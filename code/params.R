# Preset parameters for utils.R, new_server_setup.R, dada2_workflow_its.R, dada2_workflow_16s.R, and reorg_sequence_data.R

## PARAMETERS FOR NEW_SERVER_SETUP AND OUTPUT DIRECTORIES

# Parameters to define subset of data to download
# The intersection (rather than union) of these subsets is downloaded
PRESET_SITES = "all" # Can also be character vector of 4-letter NEON site codes e.g. c('ONAQ','RMNP')
PRESET_START_YR_MO = "2011-06"
PRESET_END_YR_MO = format(Sys.Date(), "%Y-%m")
PRESET_CHECK_FILE_SIZE = FALSE
PRESET_RETURN_DATA = TRUE
SAMPLE_SUBSET_PARAMS_FILENAME = "sample_subset_params.txt" # name to give file where subset params are recorded
TARGET_GENE = "ITS" # Must be "ITS", "16S", or "all"
SEQUENCING_RUNS = "all" # Can also be character vector of 5-letter NEON sequencing run IDs e.g. c('B69RF','BTP4N')
if(!(TARGET_GENE %in% c("ITS", "16S", "all"))) warning("TARGET_GENE must be 'ITS', '16S', or 'all'")

# Parameters for output directories and filenames
PRESET_OUTDIR_SEQUENCE = "/data/ZHULAB/NEON_DOB/Illumina/NEON"
PRESET_OUTDIR_SEQMETA = "/data/ZHULAB/NEON_DOB/sequence_metadata"
PRESET_OUTDIR_SOIL = "/data/ZHULAB/NEON_DOB/soil"
PRESET_OUTDIR_DADA2 = "/raid/users/claraqin/zhulab/NEON_soil_microbe_processing/data"
PRESET_FILENAME_JOINED_SEQTAB = "NEON_ITS_seqtab_nochim_DL08-13-2019.Rds"
PRESET_FILENAME_TAXTAB = "NEON_ITS_taxa_DL08-13-2019.Rds"

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

# Whether to generate additional output in dada2_workflow_its.R
VERBOSE = FALSE

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
            # (by virtue of having shorter reads)
MIN_LEN = 50 # min. allowable length of reads that pass filter
             # NOTE: it may be desirable to set this higher to increase the
             # likelihood of sufficient overlap between read pairs. However,
             # this is at odds with the incentive for setting TRUNC_Q higher
             # (see previous)
