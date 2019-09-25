# Analysis of NEON and DoB fungal sequence data
### NSF Macrosystems: Macroecology of microorganisms: Scaling fungal biodiversity from soil cores to the North American continent 

```
├── data  # NONE OF THIS DIRECTORY IS PUSHED TO GITHUB - ACCESS ON SERVER (see below)
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
└── code
    ├── data2_tutorial.R  # follows https://benjjneb.github.io/dada2/ITS_workflow.html
    |                     #  to create the processed reads from raw fastqs    
    |
    ├── downloadRawSequenceData.R  # R script for downloading NEON fastq files in bulk;
    |                              #   used to populate /data/Illumina
    |
    └── downloadAllRawSequenceData.R  # Wrapper for downloadRawSequenceData.R;
                                      #   downloads all files uploaded since given date
    
```


## Setup

**None of the data is pushed to this repo. Even if you clone this repo, you will still not be able to access the data.**

To work with the sequence data, you must create an account on `socs-stats.ucsc.edu` where it is stored. Communicate with Kai Zhu for an account. Then, do the following:

1. Log in to the server by running the Unix command `ssh [username]@socs-stats.ucsc.edu`, e.g. `ssh claraqin@socs-stats.ucsc.edu`. Enter your password when prompted.
2. Once in your home directory on the socs-stats server, run `git clone https://github.com/claraqin/NEON_DoB_analysis.git`.
3. `cd` into the new repo, e.g. `cd NEON_DoB_analysis`.
3. Establish a *symbolic link* from the central repository for sequence data (`/data/ZHULAB/NEON_DOB`) to the newly cloned repo using the `ln` Unix command: `ln -s /data/ZHULAB/NEON_DOB ./data`.

Now you should be able to work with the data on the server. There are a few ways to interact with the data. Probably the most straightforward is to use the [RStudio Server](https://socs-stats.ucsc.edu:8787). (You'll be prompted for your username and password.) You can also see [this page](https://socs-stats.ucsc.edu/doku.php) for other access options.

## Downloading raw sequence files as a batch (backend process)

To download all fastq files onto the socs-stats server, log in to the server, start a new `screen`, and run the following:

```
Rscript downloadAllRawSequenceData.R
```

(To avoid wasting time downloading files that already exist, you may have to change `startYrMo` in `downloadAllRawSequenceData.R`.)

To extract all downloaded zipped files:

```
cd /data/ZHULAB/NEON_DOB/Illumina
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

**UPDATE:** There is now a built-in function in the `neonUtilities` R package called `zipsByURI` which could replace some of this custom code.

## Pre-processing

I had to rename a single file due to a capitalization error, though there may be more in the future:

```
cd ITS
rename PLate Plate *.fastq
```

It seems that `cutadapt` (and perhaps also `filterAndTrim`) requires fastq files to be compressed.

```
gzip *.fastq
```

