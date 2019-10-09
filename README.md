# Analysis of NEON and DoB fungal sequence data
### NSF Macrosystems: Macroecology of microorganisms: Scaling fungal biodiversity from soil cores to the North American continent 

```
├── data  # NONE OF THIS DIRECTORY IS PUSHED TO GITHUB - ACCESS ON SERVER (see below)
|   |
|   ├── metadata
|   |   └── NEON_soilRawDataFiles.csv  # URL lookup table for NEON fastq downloads
|   |
|   ├── tax_ref  # taxonomic reference tables to match sequences with taxonomy
|   |   └── sh_general_release_dynamic_02.02.2019.fasta
|   |
|   └── Illumina  # raw fastq files from Illumina sequencing
|       ├── DoB
|       |   └── ITS
|       |       ├── Run1  # contains raw fastq files
|       |       ├── Run2  # contains raw fastq files
|       |       └── Run3  # contains raw fastq files
|       |   
|       └── NEON
|           ├── 16S  # contains raw fastq files
|           └── ITS  # contains raw fastq files, and directories with processed files:
|               ├── filtN  # after filtering out reads containing ambiguous bases ("N")
|               └── cutadapt2  # after removing primers/adapters
|                   └── filtered  # after passing quality filter
|                 
|
└── code  # if running R scripts in RStudio, set working directory
    |     # to be the git root directory (e.g. "NEON_DoB_analysis"), 
    |     # not "code" subdirectory
    |
    ├── utils.R  # contains various functions including one which downloads all
    |            # NEON raw microbial sequence data
    |
    ├── dada2_workflow.R  # follows https://benjjneb.github.io/dada2/ITS_workflow.html
    |                     # to process NEON ITS raw sequence data
    |
    └── dada2_to_phyloseq.R  # assembles outputs of dada2_workflow.R, plus soil data
                             # and sequence metadata, to create phyloseq object
    
```


## Setup for new collaborators

**None of the data is pushed to this repo. Even if you clone this repo, you will still not be able to access the data.**

To work with the sequence data, you must create an account on `socs-stats.ucsc.edu` where it is stored. Communicate with Kai Zhu for an account. Then, do the following:

1. Log in to the server by running the Unix command `ssh [username]@socs-stats.ucsc.edu`, e.g. `ssh claraqin@socs-stats.ucsc.edu`. Enter your password when prompted.
2. Once in your home directory on the socs-stats server, run `git clone https://github.com/claraqin/NEON_DoB_analysis.git`.
3. `cd` into the new repo, e.g. `cd NEON_DoB_analysis`.
3. Establish a *symbolic link* from the central repository for sequence data (`/data/ZHULAB/NEON_DOB`) to the newly cloned repo using the `ln` Unix command: `ln -s /data/ZHULAB/NEON_DOB ./data`.

Now you should be able to work with the data on the server. There are a few ways to interact with the data. Probably the most straightforward is to use the [RStudio Server](https://socs-stats.ucsc.edu:8787). (You'll be prompted for your username and password.) You can also see [this page](https://socs-stats.ucsc.edu/doku.php) for other access options.

## Processed data

The complete workflow thus far for processing the raw sequence data consisted of running the following scripts in order. If you need to re-run the scripts for any reason, you should set the working directory to be the git root directory (e.g. `NEON_DoB_analysis`), *not* the `code` subdirectory.

1. `dada2_workflow.R`
2. `dada2_to_phyloseq.R`

This workflow has produced the following objects of potential interest for downstream analysis:

1. a phyloseq object representing all NEON ITS data as of Aug. 13, 2019, plus their associated sample data and taxonomic table: `./code/NEON_ITS_phyloseq_DL08-13-2019.Rds` (NOTE: this phyloseq object is currently still missing soil data for over half of the samples. NEON is working to resolve this.)




To read these into `R`, use the `readRDS()` function, e.g.:

```
ps <- readRDS("./code/NEON_ITS_phyloseq_DL08-13-2019.Rds")
```


## Setup for starting from scratch on the socs-stats server

(This is how I set everything up on socs-stats in the first place.)


### Downloading raw sequence files as a batch

To download all fastq files onto the socs-stats server, log in to the server, `cd` to the root directory of your cloned repository (e.g. `cd NEON_DoB_analysis`), start a new `screen`, and start `R`. Then, in `R`:

```
source("./code/utils.R")
downloadAllRawSequenceData()
```

Notes: 
* This process may take a while. You can detach the screen while waiting, using `Ctrl+A, D`. Later you can reattach by typing `screen -r`.
* To avoid wasting time downloading files that already exist, you may have to change the `startYrMo` argument in `downloadAllRawSequenceData()`, e.g. `downloadAllRawSequenceData(startYrMo = "2019-01")`.

Once finished, quit R.

```
q()
```

To extract all downloaded zipped files:

```
cd /data/ZHULAB/NEON_DOB/Illumina/NEON
ls *.gz | xargs -n1 tar xvzf
```

To flatten the nested directories that begin with `hpc`:

```
find /data/ZHULAB/NEON_DOB/Illumina -mindepth 2 -type f -exec mv -t /data/ZHULAB/NEON_DOB/Illumina -n '{}' +
```

Then to organize ITS and 16S sequences into separate directories:

```
mkdir ITS
mv *ITS*.fastq ITS
mkdir 16S
mv *16S*.fastq 16S
```

To remove all zipped files (to save space):

```
rm *.gz
```


### Pre-processing

I had to rename a single file due to a capitalization error, though there may be more in the future:

```
cd ITS
rename PLate Plate *.fastq
```

It seems that `cutadapt` (and perhaps also `filterAndTrim`) requires fastq files to be compressed.

```
gzip *.fastq
```
