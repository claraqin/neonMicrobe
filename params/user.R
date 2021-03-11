# User parameters for the neonMicrobe package

# These parameters are used as the default values for many function arguments in the
# neonMicrobe package. Thus, they can be changed on-the-fly, but you may find it
# more convenient to set them here first.

#################################################
# GENERAL PARAMETERS

# Users can modify these parameters to suit their processing needs.

# Whether to print more detailed messages to output. This is used as the default value
# for many functions in this package, unless explicitly specified otherwise.
VERBOSE = FALSE

# Whether to use multithreading (via mclapply) in dada2 functions that allow it.
# If integer is provided instead of TRUE/FALSE, how many threads to use.
# NOTE: Windows users must set this to FALSE, because multithreading in mclapply
# is not supported on Windows systems.
MULTITHREAD = 8

#################################################
# PARAMETERS FOR SEQUENCE TRIMMING AND FILTERING

# Users can modify these parameters to suit their processing needs.

# Primers
PRIMER_ITS_FWD = "CTTGGTCATTTAGAGGAAGTAA" # Forward primer sequence
PRIMER_ITS_REV = "GCTGCGTTCTTCATCGATGC" # Reverse primer sequence
PRIMER_16S_FWD = "CCTACGGGNBGCASCAG" # Forward primer sequence
PRIMER_16S_REV = "GACTACNVGGGTATCTAATCC" # Reverse primer sequence

# Path to Cutadapt for ITS primer trimming, if necessary
CUTADAPT_PATH <- "/afs/cats.ucsc.edu/users/b/claraqin/.local/bin/cutadapt"

# Quality filtering parameters (used as defaults in dada2::filterAndTrim)
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

#################################################
# PARAMETERS FOR TAXONOMY ASSIGNMENT

# ITS: UNITE reference database (FASTA file) path
UNITE_REF_PATH = file.path(PRESET_OUTDIR_TAXREF, "sh_general_release_dynamic_04.02.2020.fasta")

# 16S: SILVA reference database (FASTA file) path
SILVA_REF_PATH = file.path(PRESET_OUTDIR_TAXREF, "silva_training_set.fa.gz")
