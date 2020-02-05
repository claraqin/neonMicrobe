# Preset parameters for utils.R, new_server_setup.R, dada2_workflow_its.R, dada2_workflow_16s.R, and reorg_sequence_data.R
preset_sites = "all"
preset_startYrMo = "2011-06"
preset_endYrMo = format(Sys.Date(), "%Y-%m")
preset_outdir_sequence = "/data/ZHULAB/NEON_DOB/Illumina/NEON"
preset_outdir_seqmeta = "/data/ZHULAB/NEON_DOB/sequence_metadata"
preset_outdir_soil = "/data/ZHULAB/NEON_DOB/soil"
preset_checkFileSize = FALSE
preset_return_data = TRUE
site_and_date_range_filename = "sites_and_date_range.txt"
# Primers
PRIMER_ITS_FWD <- "CTTGGTCATTTAGAGGAAGTAA" # Forward primer sequence
PRIMER_ITS_REV <- "GCTGCGTTCTTCATCGATGC" # Reverse primer sequence
# Cutadapt path
CUTADAPT_PATH <- "/afs/cats.ucsc.edu/users/b/claraqin/.local/bin/cutadapt"
# Whether to download and process only a small subset
SMALL_SUBSET <- TRUE # If TRUE, all fastq.gz.tar files matching these parameters will still be downloaded,
                     # and all fastq files within them will be unzipped, but only the first 2 forward-reverse
                     # pairs of sequence files from the first sequencing run ID will be processed in 
                     # dada2_workflow.R and in following scripts.