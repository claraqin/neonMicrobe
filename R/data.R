
#' 16S Sequence Metadata for NEON Great Plains Sites
#'
#' A subset of the NEON 16S sequence metadata for three
#' terrestrial NEON sites in the Great Plains: KONZ,
#' CPER, and NOGP. This sequence metadata contains
#' all records of 16S soil marker gene sequence data
#' collected from these sites in July 2017. This
#' sequence metadata has not been QCd using
#' \code{\link{qcMetadata}}.
#'
#' @format A data frame with 172 rows and 86 variables.
#' For variable definitions, see NEON soil microbe
#' marker gene sequences data product:
#' \url{https://data.neonscience.org/data-products/DP1.10108.001}
#'
#' @source
#' downloadSequenceMetadata(startYrMo = "2017-07", endYrMo = "2017-07",
#' sites = c("KONZ", "CPER", "NOGP"),
#' targetGene = "16S", outDir=FALSE)
"seqmeta_greatplains"


#' 16S Sequence Abundance Table for NEON Great Plains Sites
#'
#' A sequence abundance table generated from a
#' subset of the NEON 16S sequence metadata for three
#' terrestrial NEON sites in the Great Plains (KONZ,
#' CPER, and NOGP) collected during July 2017. This
#' sequence table has been collapsed using
#' \code{\link[dada2]{collapseNoMismatch}}.
#'
#' @details
#' This sequence abundance table was processed using maxEE=8
#' and truncLen=220 for forward and reverse reads. See the
#' full code used to generate this in the vignette
#' "Analyze NEON Great Plains 16S".
#'
#' @format A numeric matrix with 86 rows (samples) and
#' 46,730 columns (ASVs). Values represent read abundances.
#' The table has been collapsed (100% clustering) using
#' \code{\link[dada2]{collapseNoMismatch}}. Read abundances
#' are not standardized.
#'
#' @source See vignette "Analyze NEON Great Plains 16S".
"seqtab_greatplains"


#' Soil Physical and Chemical Properties for NEON Great Plains Sites
#'
#' A data frame containing soil physical and chemical
#' properties data (DP1.10086.001: "Soil physical and chemical properties,
#' periodic", tables sls_soilCoreCollection, sls_soilMoisture, sls_soilpH,
#' and sls_soilChemistry) for three terrestrial NEON sites in the Great
#' Plains (KONZ, CPER, and NOGP) collected during July 2017.
#'
#' @details
#' This sequence abundance table was processed using maxEE=8
#' and truncLen=220 for forward and reverse reads. See the
#' full code used to generate this in the vignette
#' "Analyze NEON Great Plains 16S".
#'
#' @format A numeric matrix with 86 rows (samples) and
#' 46,730 columns (ASVs). Values represent read abundances.
#' The table has been collapsed (100% clustering) using
#' \code{\link[dada2]{collapseNoMismatch}}. Read abundances
#' are not standardized.
#'
#' @source See vignette "Analyze NEON Great Plains 16S".
"soil_greatplains"
