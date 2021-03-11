# System parameters for the neonMicrobe package

#################################################
# PARAMETERS FOR FILE STRUCTURE SETUP

# WE GENERALLY RECOMMEND AGAINST MODIFYING THESE PARAMETERS EXCEPT FOR CODE DEVELOPMENT PURPOSES.

# Main directory structure for the pipeline. Creates nested directory structure for inputs and outputs.
# Can leave the following as-is to generate default directory structure.
# Or, character string to generate and use alternative directory structure.
PRESET_OUTDIR = file.path(getwd(), "data/")
PRESET_OUTDIR_SEQUENCE = file.path(PRESET_OUTDIR, "raw_sequence") # for sequence data (fastq files)
PRESET_OUTDIR_SEQMETA = file.path(PRESET_OUTDIR, "sequence_metadata") # for sequence metadata
PRESET_OUTDIR_SOIL = file.path(PRESET_OUTDIR, "soil") # for soil data
PRESET_OUTDIR_TAXREF = file.path(PRESET_OUTDIR, "tax_ref") # for taxonomy reference data
PRESET_OUTDIR_OUTPUTS = file.path(getwd(), "outputs") # for outputs (sequence table, taxonomy table, phyloseq object)
