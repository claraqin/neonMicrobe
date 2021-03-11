# Parameter validity tests for the neonMicrobe package

#################################################
# VALIDITY CHECKS

## DO NOT MODIFY EXCEPT FOR CODE DEVELOPMENT PURPOSES.

# If Windows users attempt to set MULTITHREAD to anything but FALSE, print a warning
if (Sys.info()["sysname"] == "Windows"){
  if(!identical(MULTITHREAD, FALSE)){
    warning("On Windows systems, multithreading (via mclapply) is not available. Setting MULTITHREAD to FALSE.")
    MULTITHREAD = FALSE
  }
}

# If CUTADAPT_PATH is invalid, print a warning
if(!file.exists(CUTADAPT_PATH)){
  warning("Could not find the Cutadapt tool in file system. This tool is necessary for ITS (fungal) raw sequence processing. Please load the module if necessary on an HPC (i.e. 'module load cutadapt') or download from https://cutadapt.readthedocs.io/en/stable/installation.html")
}

# If UNITE_REF_PATH is invalid, print a warning
if(!file.exists(UNITE_REF_PATH)){
  warning("Could not find the UNITE reference database at the provided filepath: ", UNITE_REF_PATH)
}

# If SILVA_REF_PATH is invalid, print a warning
if(!file.exists(SILVA_REF_PATH)){
  warning("Could not find the SILVA reference database at the provided filepath: ", SILVA_REF_PATH)
}


# if(!(TARGET_GENE %in% c("ITS", "16S", "all"))) warning("TARGET_GENE must be 'ITS', '16S', or 'all'")
