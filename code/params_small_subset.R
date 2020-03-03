# Preset parameters for utils.R, new_server_setup.R, dada2_workflow_its.R, dada2_workflow_16s.R, and reorg_sequence_data.R
preset_site = "ONAQ"
preset_startYrMo = "2017-10"
preset_endYrMo = "2017-10"
preset_outdir_sequence = "/data/ZHULAB/NEON_DOB/Illumina/NEON/small_subset"
preset_outdir_seqmeta = "/data/ZHULAB/NEON_DOB/sequence_metadata/small_subset"
preset_outdir_soil = "/data/ZHULAB/NEON_DOB/soil/small_subset"
preset_checkFileSize = FALSE
preset_return_data = TRUE
site_and_date_range_filename = "sites_and_date_range.txt"
# Primers
PRIMER_ITS_FWD <- "CTTGGTCATTTAGAGGAAGTAA" # Forward primer sequence
PRIMER_ITS_REV <- "GCTGCGTTCTTCATCGATGC" # Reverse primer sequence
# Path to the cutadapt command
CUTADAPT_PATH <- "/afs/cats.ucsc.edu/users/b/claraqin/.local/bin/cutadapt"
# Whether to process only a small subset
SMALL_SUBSET <- TRUE # If TRUE, all fastq.gz.tar files matching these parameters will still be downloaded,
                     # and all fastq files within them will be unzipped, but only the first 2 forward-reverse
                     # pairs of sequence files from the first sequencing run ID will be processed in 
                     # dada2_workflow_its.R and in following scripts.