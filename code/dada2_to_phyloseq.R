# Merge various NEON data products to link the outputs of the
# DADA2 pipeline (e.g. non-chimeric sequence table, taxonomic table)
# to environmental/sample data.

source("../../code/utils.R")

metadata <- downloadAllSequenceMetadata()

metadata$geneticSampleID <- sub("-DNA[1-3]", "", metadata$dnaSampleID)

# TODO: Merge metadata with the output of downloadAllRawSoilData(),
# joining by geneticSampleID. Then match rownames of seqtab_nochim
# to the environmental data via the newly merged dataset's 
# "internalLabID" column.

