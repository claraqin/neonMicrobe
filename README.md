# Analysis of NEON and DoB fungal sequence data
### NSF Macrosystems: Macroecology of microorganisms: Scaling fungal biodiversity from soil cores to the North American continent 

```
├── data
|   ├── metadata
|   |   └── NEON_soilRawDataFiles.csv  # URL lookup table for NEON fastq downloads
|   |
|   └── Illumina  # raw fastq files from NEON and DoB Illumina sequencing;
|                 #   NOT PUSHED TO GITHUB - ACCESS ON SERVER (see below)
|
└── code
    ├── downloadRawSequenceData.R  # R script for downloading NEON fastq files in bulk;
    |                              #   used to populate /data/Illumina
    |
    └── downloadAllRawSequenceData.R  # Wrapper for downloadRawSequenceData.R;
                                      #   downloads all files uploaded since given date
    
```


## Setup

**None of the raw sequence data is pushed to this repo. If you clone this repo, you will not be able to access the raw sequence data.**

To work with the sequence data, you must create an account on `socs-stats.ucsc.edu` where it is stored. Communicate with Kai for an account. Then, do the following:

1. Log in to the server by running the Unix command `ssh [username]@socs-stats.ucsc.edu`, e.g. `ssh claraqin@socs-stats.ucsc.edu`. Enter your password when prompted.
2. Once in your home directory on the socs-stats server, run `git clone https://github.com/claraqin/NEON_DoB_analysis.git`.
3. `cd` into the new repo, e.g. `cd NEON_DoB_analysis`.
3. Establish a *symbolic link* from the central repository for sequence data (`/data/ZHULAB/NEON_DOB/Illumina`) to the newly cloned repo using the `ln` Unix command: `ln -s /data/ZHULAB/NEON_DOB/Illumina ./data/Illumina`.

Now you should be able to work with the raw sequence data on the server. There are a few ways to interact with the data. Probably the most straightforward is to use the [RStudio Server](https://socs-stats.ucsc.edu:8787). (You'll be prompted for your username and password.) You can also see [this page](https://socs-stats.ucsc.edu/doku.php) for other access options.

## Downloading raw sequence files as a batch

To download all fastq files, run on `screen`:

```
Rscript downloadAllRawSequenceData.R
```

(To avoid wasting time downloading files that already exist, you may have to change `startYrMo` in `downloadAllRawSequenceData.R`.)

To extract all downloaded zipped files:

```
cd /data/ZHULAB/NEON_DOB/Illumina
ls *.gz | xargs -n1 tar xvzf
```