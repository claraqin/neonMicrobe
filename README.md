# neonMicrobe

`neonMicrobe` is a suite of functions for downloading, pre-processing, and assembling heterogeneous data around the NEON soil microbe marker gene sequence data. To do so, `neonMicrobe` downloads NEON data products from the NEON Data API and processes sequences using the DADA2 workflow. In addition, `neonMicrobe` offers a processing-batch infrastructure to encourage explicit versioning of processed data.

## Installation

The development version of `neonMicrobe` can be installed directly from this GitHub repo using this code:

```
install.packages("devtools")
devtools::install_github("claraqin/neonMicrobe")
```

## User-installed dependencies

In addition to the R package dependencies which are installed alongside `neonMicrobe`, users may also need to complete the following requirements before using some functions in `neonMicrobe`:

1. **For taxonomic assignment via DADA2**, you will need to install the latest taxonomic reference datasets for ITS or 16S sequences. Consult the [DADA2 taxonomic reference data webpage](https://benjjneb.github.io/dada2/training.html) for more information. For organizational purposes, we recommend keeping these files in the `data/tax_ref` subdirectory that is created after you run `makeDataDirectories()` (see "Input data" below).
2. **For trimming of ITS sequence primers**, you will need to install `cutadapt`. Installation instructions can be found [here](https://cutadapt.readthedocs.io/en/stable/installation.html). Once installed, you can tell `neonMicrobe` where to look for it by specifying the `cutadapt_path` argument each time you use the `trimPrimerITS` function. For an example, see the "Process 16S Sequences" vignette or the "Process ITS Sequences" vignette.

## Overview

Tutorials for `neonMicrobe` are available in the `vignettes` directory.

1. **Download NEON Data** – how to use the functions in this package to download the specific scope NEON soil microbe marker gene sequence data and associated data relevant to your analysis. Leverages the `neonUtilities` [R package](https://github.com/NEONScience/NEON-utilities).
2. **Process Sequences** – how to use the functions in this package to add associated environmental variables to the ASV tables. Leverages the `dada2` [R package](https://github.com/benjjneb/dada2). The dada denoising algorithm partitions reads into amplicon sequence variants (ASVs), which are finer in resolution than OTUs.
3. **Add Environmental Variables** – how to use the functions in this package to add associated environmental variables to the ASV tables. Joins the data together in the form of one or more Phyloseq objects.
4. (Optional) **Sensitivity Analysis** – how to use the functions in this package to test the effects of quality filtering parameters and decisions on the resulting ecological inference. Can be used to test for an acceptable range of custom parameters.

![NEON Ecosphere MS Figure-Making Workspace (3)](https://user-images.githubusercontent.com/12421420/111393342-ce937d00-8675-11eb-8b63-530aced18352.png)

## Data storage structure

### Input data

The **Download NEON Data** vignette demonstrates how to download NEON data, optionally writing to the file system. By default, the input data is downloaded into the following structure, which is created in the working directory after running `makeDataDirectories()`:

![NEON Ecosphere MS Figure-Making Workspace (10)](https://user-images.githubusercontent.com/12421420/113089173-f3badc00-919b-11eb-84e6-b7f9a2abbb72.png)

The tree structure in the upper-left represents the data directory structure constructed within the project root directory. Red dotted lines represent explicit linkages between NEON data products via shared data fields. (a) Sequence metadata is downloaded from NEON data product DP1.10108.001 (Soil microbe marker gene sequences) using the downloadSequenceMetadata() function. (b) Raw microbe marker gene sequence data is downloaded from NEON based on the sequence metadata using the downloadRawSequenceData() function. (c) Soil physical and chemical data is downloaded from NEON data product DP1.10086.001 using the downloadRawSoilData() function. (d) Taxonomic reference datasets (e.g. SILVA, UNITE) are added separately by the user.

### Output data

The **Process Sequences** and **Add Environmental Variables** vignettes demonstrate how to process the NEON data inputs into useful sample-abundance tables with accompanying environmental data.

(Under construction.)


## Methods paper

A methods paper describing the use of `neonMicrobe` is currently in review.

## Examples

### Quick start

[Analyze NEON Great Plains 16S Data](https://people.ucsc.edu/~claraqin/analyze-neon-greatplains-16s.R)

### Vignettes

1. [Download NEON Data](https://people.ucsc.edu/~claraqin/download-neon-data.html)
2. [Process 16S Sequences](https://people.ucsc.edu/~claraqin/process-16s-sequences.html) or [Process ITS Sequences](https://people.ucsc.edu/~claraqin/process-its-sequences.html)
3. [Add Environmental Variables to 16S Data](https://people.ucsc.edu/~claraqin/add-environmental-variables-16s.html)
