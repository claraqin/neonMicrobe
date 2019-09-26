# Download ALL raw soil data
#
# Data product are
# - DP1.10078: "Soil chemical properties (Distributed periodic)"
# - DP1.10086: "Soil physical properties (Distributed periodic)"
#
# Wrapper for downloadRawSoilData.R,
# parameterized for socs-stats.ucsc.edu server,
# and will pull data up to current month.

library(neonUtilities)

source("downloadRawSoilData.R")

sites = "all"
startYrMo = "2011-06"
endYrMo = format(Sys.Date(), "%Y-%m")
outdir = "/afs/cats.ucsc.edu/users/b/claraqin/zhulab/NEON_DoB_analysis/data/metadata"
checkFileSize = FALSE
return_data = TRUE

soil <- downloadRawSoilData(sites, startYrMo, endYrMo, outdir, checkFileSize, return_data)
