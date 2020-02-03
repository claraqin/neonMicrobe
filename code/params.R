# Preset parameters for the "All" variants of download functions:
# (Also used in dada2_workflow_its.R, dada2_workflow_16s.R, and reorg_sequence_data.R)
preset_sites = "all"
preset_startYrMo = "2011-06"
preset_endYrMo = format(Sys.Date(), "%Y-%m")
preset_outdir_sequence = "/data/ZHULAB/NEON_DOB/Illumina/NEON/raw3"
preset_outdir_seqmeta = "/data/ZHULAB/NEON_DOB/sequence_metadata/test"
preset_outdir_soil = "/data/ZHULAB/NEON_DOB/soil"
preset_checkFileSize = FALSE
preset_return_data = TRUE
site_and_date_range_filename = "sites_and_date_range.txt"
# Primers
PRIMER_ITS_FWD <- "CTTGGTCATTTAGAGGAAGTAA" # Forward primer sequence
PRIMER_ITS_REV <- "GCTGCGTTCTTCATCGATGC" # Reverse primer sequence
# Cutadapt path
CUTADAPT_PATH <- "/afs/cats.ucsc.edu/users/b/claraqin/.local/bin/cutadapt"