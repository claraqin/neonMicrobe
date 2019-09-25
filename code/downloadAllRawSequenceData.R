# Download ALL raw sequence data
#
# Metadata is DP1.10108: "Soil microbe marker gene sequences"
#
# Wrapper for downloadRawSequenceData.R,
# parameterized for socs-stats.ucsc.edu server,
# and will pull data up to current month.

library(neonUtilities)
library(utils)

source("downloadRawSequenceData.R")

sites = "all"
startYrMo = "2011-06"
endYrMo = format(Sys.Date(), "%Y-%m")
outdir = "/afs/cats.ucsc.edu/users/b/claraqin/zhulab/NEON_DoB_analysis/data/Illumina"
checkFileSize = FALSE
return_data = FALSE

downloadRawSequenceData(sites, startYrMo, endYrMo, outdir, checkFileSize,
                        return_data)
