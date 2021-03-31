# neonMicrobe

`neonMicrobe` is a suite of functions for downloading, pre-processing, and assembling heterogeneous data around the NEON soil microbe marker gene sequence data. To do so, `neonMicrobe` downloads NEON data products from the NEON Data API and processes sequences using the DADA2 workflow. In addition, `neonMicrobe` offers a processing-batch infrastructure to encourage explicit versioning of processed data.

## Installation

`neonMicrobe` is currently still in development, but it will soon be installable from GitHub using this code:

```
install.packages("devtools")
devtools::install_github("claraqin/neonMicrobe")
```

## User-installed dependencies

In addition to the R package dependencies which are installed alongside `neonMicrobe`, users may also need to complete the following requirements before using some functions in `neonMicrobe`:

1. **For taxonomic assignment via DADA2**, you will need to install the latest taxonomic reference datasets for ITS or 16S sequences. Consult the [DADA2 taxonomic reference data webpage](https://benjjneb.github.io/dada2/training.html) for more information. The `params/user.R` file contains optional parameter slots for specifying the file path to the UNITE and SILVA reference datasets, respectively.
2. **For trimming of ITS sequence primers**, you will need to install `cutadapt`. Installation instructions can be found [here](https://cutadapt.readthedocs.io/en/stable/installation.html). Once installed, you can tell `neonMicrobe` where to look for it by either:
    a. changing the CUTADAPT_PATH variable in `params/user.R` to match the location of the `cutadapt` file on your system, or
    b. specifying the `cutadapt_path` argument each time you use the `trimPrimerITS` function.

## Overview

Tutorials for `neonMicrobe` are available in the `vignettes` directory.

1. **Download NEON Data** – how to use the functions in this package to download the specific scope NEON soil microbe marker gene sequence data and associated data relevant to your analysis. Leverages the `neonUtilities` [R package](https://github.com/NEONScience/NEON-utilities).
2. **Process Sequences** – how to use the functions in this package to add associated environmental variables to the ASV tables. Leverages the `dada2` [R package](https://github.com/benjjneb/dada2). The dada denoising algorithm partitions reads into amplicon sequence variants (ASVs), which are finer in resolution than OTUs.
3. **Add Environmental Variables** – how to use the functions in this package to add associated environmental variables to the ASV tables. Joins the data together in the form of one or more Phyloseq objects.
4. (Optional) **Sensitivity Analysis** – how to use the functions in this package to test the effects of quality filtering parameters and decisions on the resulting ecological inference. Can be used to test for an acceptable range of custom parameters.

![NEON Ecosphere MS Figure-Making Workspace (3)](https://user-images.githubusercontent.com/12421420/111393342-ce937d00-8675-11eb-8b63-530aced18352.png)

## User-specified parameters

Users should feel empowered to specify custom parameters in the `params/user.R`. The functions in `neonMicrobe` will refer to these parameters to set default values for many function arguments. Of course, these arguments can be changed on-the-fly, but you may find it more convenient to set them in `params/user.R` first. Parameters include but may not be limited to the following:

- VERBOSE
- MULTITHREAD
- PRIMER_ITS_FWD
- PRIMER_ITS_REV
- PRIMER_16S_FWD
- PRIMER_16S_REV
- CUTADAPT PATH
- MAX_EE_FWD
- MAX_EE_REV
- TRUNC_Q
- MIN_LEN
- UNITE_REF_PATH
- SILVA_REF_PATH

See `params/user.R` for more details.


## Data storage structure

### Input data

The **Download NEON Data** vignette demonstrates how to download NEON data, optionally writing to the file system. By default, the input data is downloaded into the following structure:

![NEON Ecosphere MS Figure-Making Workspace (4)](https://user-images.githubusercontent.com/12421420/111392531-05689380-8674-11eb-80ee-92bb009889c9.png)

The tree structure in the upper-left represents the data directory structure constructed within the project root directory. Red dotted lines represent explicit linkages between NEON data products via shared data fields. (a) Sequence metadata is downloaded from NEON data product DP1.10108.001 (Soil microbe marker gene sequences) using the downloadSequenceMetadata() function. (b) Raw microbe marker gene sequence data is downloaded from NEON based on the sequence metadata using the downloadRawSequenceData() function. (c) Soil physical and chemical data is downloaded from NEON data product DP1.10086.001 using the downloadRawSoilData() function. (d) Taxonomic reference datasets (e.g. SILVA, UNITE) are added separately by the user.

### Output data

The **Process Sequences** and **Add Environmental Variables** vignettes demonstrate how to process the NEON data inputs into useful sample-abundance tables with accompanying environmental data.

(Under construction.)


## Methods paper

A methods paper describing the use of `neonMicrobe` is currently in preparation.
