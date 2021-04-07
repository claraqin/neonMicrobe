# neonMicrobe

`neonMicrobe` is a suite of functions for downloading, pre-processing, and assembling heterogeneous data around the NEON soil microbe marker gene sequence data. To do so, `neonMicrobe` downloads NEON data products from the NEON Data API and processes sequences using the DADA2 workflow. In addition, `neonMicrobe` offers a processing-batch infrastructure to encourage explicit versioning of processed data.

## Installation

**Temporary solution for beta-testers until this repo is made public:**

1. Download the package bundle [here](https://drive.google.com/file/d/1jKvkSdbu23gbXMDMjbiaARHEHG74iV12/view?usp=sharing). It should be `neonMicrobe_0.0.0.9000.tar.gz`. Drag it into a directory that you can easily navigate to in R.
2. In R, run the following:

```
setwd([directory containing the package bundle])
install.packages("neonMicrobe_0.0.0.9000.tar.gz", repos = NULL, type = "source")
```

The development version of `neonMicrobe` can be installed directly from this GitHub repo using this code:

```
install.packages("devtools")
devtools::install_github("claraqin/neonMicrobe")
```

## User-installed dependencies

In addition to the R package dependencies which are installed alongside `neonMicrobe`, users may also need to complete the following requirements before using some functions in `neonMicrobe`:

1. **For taxonomic assignment via DADA2**, you will need to install the latest taxonomic reference datasets for ITS or 16S sequences. Consult the [DADA2 taxonomic reference data webpage](https://benjjneb.github.io/dada2/training.html) for more information. For organizational purposes, we recommend keeping these files in the `data/tax_ref` subdirectory that is created after you run `makeDataDirectories()` (see "Input data" below).
2. **For trimming of ITS sequence primers**, you will need to install `cutadapt`. Installation instructions can be found [here](https://cutadapt.readthedocs.io/en/stable/installation.html). Once installed, you can tell `neonMicrobe` where to look for it by specifying the `cutadapt_path` argument each time you use the `trimPrimerITS` function. For an example, see the "Process 16S Sequences" vignette or the "Process ITS Sequences" vignette.

## Quick start

The following R script makes use of `neonMicrobe` to create ASV tables for 16S sequences collected from three NEON sites in the Great Plains:

[Analyze NEON Great Plains 16S Sequences](https://people.ucsc.edu/~claraqin/analyze-neon-greatplains-16s.R)

## Overview

Tutorials for `neonMicrobe` are available in the `vignettes` directory, and some are also linked here:

1. **[Download NEON Data](https://people.ucsc.edu/~claraqin/download-neon-data.html)** – how to use the functions in this package to download the specific scope NEON soil microbe marker gene sequence data and associated data relevant to your analysis. Leverages the `neonUtilities` [R package](https://github.com/NEONScience/NEON-utilities).
2. **[Process 16S Sequences](https://people.ucsc.edu/~claraqin/process-16s-sequences.html)** and **[Process ITS Sequences](https://people.ucsc.edu/~claraqin/process-its-sequences.html)** – how to use the functions in this package to add associated environmental variables to the ASV tables. Leverages the `dada2` [R package](https://github.com/benjjneb/dada2). The dada denoising algorithm partitions reads into amplicon sequence variants (ASVs), which are finer in resolution than OTUs.
3. **[Add Environmental Variables to 16S Data](https://people.ucsc.edu/~claraqin/add-environmental-variables-16s.html)** – how to use the functions in this package to add associated environmental variables to the ASV tables. Joins the data together in the form of one or more Phyloseq objects.
4. (Optional) **Sensitivity Analysis** – how to use the functions in this package to test the effects of quality filtering parameters and decisions on the resulting ecological inference. Can be used to test for an acceptable range of custom parameters.
5. (Coming soon) **Processing Batches** - how to use the processing-batch feature to keep track of the parameters used to create various sets of output data.

![NEON Ecosphere MS Figure-Making Workspace (3)](https://user-images.githubusercontent.com/12421420/111393342-ce937d00-8675-11eb-8b63-530aced18352.png)

## File storage structure

### Input data

The **[Download NEON Data](https://people.ucsc.edu/~claraqin/download-neon-data.html)** vignette demonstrates how to download NEON data, optionally writing to the file system. By default, the input data is downloaded into the following structure, which is created in the working directory after running `makeDataDirectories()`:

![NEON Ecosphere MS Figure-Making Workspace (10)](https://user-images.githubusercontent.com/12421420/113089173-f3badc00-919b-11eb-84e6-b7f9a2abbb72.png)

The tree structure in the upper-left represents the data directory structure constructed within the project root directory. Red dotted lines represent explicit linkages between NEON data products via shared data fields. (a) Sequence metadata is downloaded from NEON data product DP1.10108.001 (Soil microbe marker gene sequences) using the downloadSequenceMetadata() function. (b) Raw microbe marker gene sequence data is downloaded from NEON based on the sequence metadata using the downloadRawSequenceData() function. (c) Soil physical and chemical data is downloaded from NEON data product DP1.10086.001 using the downloadRawSoilData() function. (d) Taxonomic reference datasets (e.g. SILVA, UNITE) are added separately by the user.

### Output data

The **Process ([16S](https://people.ucsc.edu/~claraqin/process-16s-sequences.html)/[ITS](https://people.ucsc.edu/~claraqin/process-its-sequences.html)) Sequences** and **[Add Environmental Variables to 16S Data](https://people.ucsc.edu/~claraqin/add-environmental-variables-16s.html)** vignettes demonstrate how to process the NEON data inputs into useful sample-abundance tables with accompanying environmental data.

By default, output data from `neonMicrobe` is downloaded into the `outputs/` directory.

```
─ outputs
  ├── mid_process
  │   ├── 16S
  │   └── ITS
  └── track_reads
      ├── 16S
      └── ITS
```

The `mid_process/` subdirectory contains files in the middle of being processed -- for example, fastq files that have been trimmed or filtered, and sequencing run-specific ASV tables that have not yet been joined together. Once the desired outputs have been created, you may choose to clear the contents of `mid_process/`, or leave them to retrace your processing steps. 

The `track_reads/` subdirectory contains tables tracking the number of reads remaining at each step in the pipeline, from the "raw" sequence files downloaded from NEON to the ASV table. These tables can be useful for pinpointing steps and samples for which an unusual number of reads were lost.

(Coming soon: When the processing batch feature is released, the default outputs directory will be switched to `batch_outputs`. More on this later!)

## Methods paper

A methods paper describing the use of `neonMicrobe` is currently in review.
