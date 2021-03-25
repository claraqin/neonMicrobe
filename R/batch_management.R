# Create new environments for storing processing batch parameters

# Package environment structure (builds on dada_opts)
#
# dada_opts
# └── neonmicrobe_env
#     ├── batch_env1
#     ├── batch_env2
#     └── batch_env3 ...


valid_params <- c()

#' Assign Final Value
#'
#' Helper function. Assigns an immutable value to a variable name within a specified
#' environment. This function is meant to be used in the context of setting processing
#' batch parameters.
#'
#' @param name Variable name to assign value to.
#' @param value Value assigned.
#' @param env Environment in which to conduct the assignment. Defaults to the calling environment.
#' @param verbose Default TRUE. Whether to print warning message if a parameter has already been assigned.
#'
#' @return No value is returned.
#'
#' @examples
assign_final <- function(name, value, env = parent.frame(), verbose=TRUE) {
  if(exists(name, envir=env, inherits=FALSE)) {
    if(verbose) {
      warning("Cannot set parameter '", name, "' because it has already been set within this processing batch. ",
              "To start a new batch, use newBatch().")
    }
  } else {
    assign(name, value, envir=env)
  }
}

#' Load Processing Parameters from File
#'
#' Loads neonMicrobe processing parameters from a specified file to a specified environment,
#' usually associated with a processing batch.
#'
#' @param params_file Filepath to R script containing name-value assignments (e.g. x = 1) that will be used to create an environment.
#' @param env Environment in which to load the parameters. By default (NULL), searches for current processing batch environment. If found, loads parameters there.
#' @param verbose Default TRUE. Whether to print warning message if a parameter has already been assigned in a processing batch.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' e <- new.env()
#' ls(e)  # character(0)
#' loadParameters("./R/test_params.R", env=e)
#' ls(e)  # "x" "y"
#' get("x", e)  # "hello there"
#' get("y", e)  # "goodbye now"
#' loadParameters("./R/test_params.R", env=e) # warning
#' }
loadParameters <- function(params_file, env = NULL, verbose=TRUE) {
  if(is.null(env)) {
    if("WORKING_BATCH_ID" %in% ls(envir = neonmicrobe_env)) {
      message("Processing batch found. Loading parameters into run ID '", get("WORKING_BATCH_ID", neonmicrobe_env), "'.")
      env <- batch_env
    } else {
      warning("Processing batch is not currently set. Use setBatch() to do so.")
      return(invisible(NULL))
    }
  }
  source(params_file, local=TRUE)
  params <- ls()
  params <- params[!(params %in% formalArgs(loadParameters))]
  for(p in params) {
    assign_final(p, get(p), env=env, verbose=verbose)
  }
}


#' Save Processing Parameters to File
#'
#' Saves neonMicrobe processing parameters from a specified environment to a specified file,
#' usually associated with a processing batch.
#'
#' @param params_file (Optional) Filepath to R script in which the environment will be saved as text of name-value assignments (e.g. x = 1). By default, uses file path in the current batch directory called "params.R".
#' @param env Environment containing parameters to be saved. By default (NULL), searches for current processing batch environment. If found, saves parameters from there.
#' @param verbose Default TRUE. Whether to print message confirming that parameters were saved.
#'
#' @return
#' @export
#'
#' @examples
saveParameters <- function(params_file = NULL, env = NULL, verbose = TRUE) {
  if (is.null(params_file)) {
    if ("WORKING_BATCH_ID" %in% ls(envir = neonmicrobe_env)) {
      params_file <- file.path(get("WORKING_BATCH_DIR", neonmicrobe_env), "params.R")
    } else {
      warning("saveParameters: Invalid params_file: Processing batch is not currently set. Use setBatch() to do so.")
      return(invisible(NULL))
    }
  }
  if(is.null(env)) {
    if ("WORKING_BATCH_ID" %in% ls(envir = neonmicrobe_env)) {
      env <- batch_env
    } else {
      warning("saveParameters: Invalid env: Processing batch is not currently set. Use setBatch() to do so.")
      return(invisible(NULL))
    }
  }
  params <- as.list.environment(env)
  for(p_nm in names(params)) {
    if(is.character(params[[p_nm]])) {
      params[[p_nm]] <- shQuote(params[[p_nm]])
    }
  }
  params_con <- file(params_file)
  if(length(params) > 0) {
    writeLines(paste(names(params), "=", params), params_con)
  } else {
    writeLines("", params_con)
  }

  close(params_con)
  if(verbose) message("Processing parameters saved to ", params_file)
}


#' Create New Processing Batch
#'
#' Creates and switches to a new processing batch for enhancing the
#' reproducibility of NEON soil microbe marker gene sequence workflows.
#' Set batch processing parameters with \code{\link{setBatchParameter}}, or
#' view them with \code{\link{getBatchParameter}}.
#'
#' @param seqmeta_file Character string. File path to the sequence metadata to associate with this processing batch. Once set, this cannot be changed except by overwriting this batch.
#' @param batch_id Character string. Unique ID to use for the new processing batch.
#' @param batches_dir Default file.path(NEONMICROBE_DIR_OUTPUTS(), "batches"). Directory where batch-specific directories are saved.
#' @param overwrite Default FALSE. If processing batch already exists in the specified directory, whether to overwrite it.
#'
#' @return No value returned
#' @export
#'
#' @examples
#' \dontrun{
#' newBatch("data/sequence_metadata/mmg_soilMetadata_ITS_2021-03-08134134.csv") # creates new batch based on timestamp
#' newBatch("data/sequence_metadata/mmg_soilMetadata_ITS_2021-03-08134134.csv", batch_id="abc") # creates new batch based on user-specified name
#' }
newBatch <- function(seqmeta_file, batch_id = NULL, batches_dir = file.path(NEONMICROBE_DIR_OUTPUTS(), "batches"), overwrite=FALSE) {
  if(!dir.exists(batches_dir)) dir.create(batches_dir, recursive=TRUE)
  if(is.numeric(batch_id)) batch_id <- as.character(batch_id)

  # If batch_id is not provided, use the current time
  if(is.null(batch_id)) {
    batch_id <- sub(" ", "_", gsub(":", "", Sys.time()))
    while(dir.exists(file.path(batches_dir, batch_id))) {
      Sys.sleep(1)
      batch_id <- sub(" ", "_", gsub(":", "", Sys.time()))
    }
  }

  # Create batch
  if(!file.exists(seqmeta_file)) {
    warning("Processing batch could not be created: No sequence metadata found at specified file location.")
  } else {
    new_batch_dir <- file.path(batches_dir, batch_id)

    if(dir.exists(new_batch_dir) & !identical(overwrite, TRUE)) {
      warning("Directory for processing batch ID ", batch_id, " already exists: ", new_batch_dir, ". To really overwrite, set overwrite=TRUE.")
    } else {

      # Overwrite if batch already exists and overwrite==TRUE
      if(dir.exists(new_batch_dir) & identical(overwrite, TRUE)) unlink(new_batch_dir, recursive=TRUE)

      # Create new batch dir
      dir.create(new_batch_dir)

      # Make ABOUT file in batch directory
      about_file <- file(file.path(new_batch_dir, "ABOUT"))
      writeLines(c(paste0("Batch created: ", Sys.time()),
                   paste0("neonMicrobe version: ", packageVersion("neonMicrobe")),
                   paste0("Sequence metadata: ", normalizePath(seqmeta_file))),
                 about_file)
      close(about_file)

      # Make directory structure in new batch directory
      createDirIfNotExist <- function(dir) {
        if(!dir.exists(dir)) dir.create(dir, recursive=TRUE)
      }

      # First, create intermediary directories for ITS and 16S data in the middle
      # of being processed
      processing_its_dir <- file.path(new_batch_dir, "mid_process", "ITS")
      processing_16s_dir <- file.path(new_batch_dir, "mid_process", "16S")
      createDirIfNotExist(file.path(processing_its_dir, "1_filtN"))
      createDirIfNotExist(file.path(processing_its_dir, "2_trimmed"))
      createDirIfNotExist(file.path(processing_its_dir, "3_filtered"))
      createDirIfNotExist(file.path(processing_its_dir, "4_seqtabs"))
      createDirIfNotExist(file.path(processing_16s_dir, "1_trimmed"))
      createDirIfNotExist(file.path(processing_16s_dir, "2_filtered"))
      createDirIfNotExist(file.path(processing_16s_dir, "3_seqtabs"))

      # Also create directories for read-tracking tables
      read_tracking_its_dir <- file.path(new_batch_dir, "track_reads", "ITS")
      read_tracking_16s_dir <- file.path(new_batch_dir, "track_reads", "16S")
      createDirIfNotExist(read_tracking_its_dir)
      createDirIfNotExist(read_tracking_16s_dir)

      message("Created new batch directory at ", new_batch_dir)

      # Switch to the new batch
      setBatch(batch_id, batches_dir, suppress_load_parameters = TRUE)

      # Set the parameter for the batch-specific outputs directory
      setBatchParameter(DIR_OUTPUTS = normalizePath(new_batch_dir), verbose=FALSE)
      setBatchParameter(SEQMETA_FILE = normalizePath(seqmeta_file), verbose=FALSE)

      # Save parameters to local params.R
      saveParameters(file.path(new_batch_dir, "params.R"), batch_env, verbose=FALSE)
    }
  }
}


#' Get Current Processing Batch ID
#'
#' Get the unique ID of the current processing batch.
#'
#' @param verbose Default FALSE. If TRUE, returns information about the current processing batch.
#'
#' @return Character. Unique ID of the current processing batch.
#' @export
#'
#' @examples
#' getBatch() # "xyz"
#' getBatch(verbose = TRUE) # "xyz", and table containing info about batch "xyz"
getBatch <- function(verbose=FALSE) {
  if ("WORKING_BATCH_ID" %in% ls(envir = neonmicrobe_env)) {
    batch_id <- get("WORKING_BATCH_ID", envir = neonmicrobe_env)
    if(verbose==TRUE) {
      this_batch_dir <- get("WORKING_BATCH_DIR", envir = neonmicrobe_env)
      about <- readLines(con = file.path(this_batch_dir, "ABOUT"), 3)
      print(about)
    }
    return(batch_id)
  } else {
    warning("Processing batch is not currently set. Use setBatch() to do so.")
    return(invisible(NULL))
  }
}

#' Switch to Existing Processing Batch
#'
#' @param batch_id Character string. Unique ID of an existing processing batch.
#' @param batches_dir Default file.path(NEONMICROBE_DIR_OUTPUTS(), "batches"). Directory where batch-specific directories are saved.
#' @param suppress_load_parameters Default FALSE. If TRUE, does not attempt to load parameters into the batch environment. Intended for use in newBatch()
#' @param verbose Default TRUE. Whether to print message confirming that batch has been set.
#'
#' @return
#' @export
#'
#' @examples
setBatch <- function(batch_id = NULL, batches_dir = file.path(NEONMICROBE_DIR_OUTPUTS(), "batches"),
                     suppress_load_parameters = FALSE, verbose = TRUE) {
  if(!is.null(batch_id) & !identical("character", class(batch_id))) {
    stop("Batch ID provided was not a character string. Unable to set batch.")
  }
  # If batch_id is NOT provided, clear the WORKING_BATCH_ID variable and delete the batch-specific env
  if(is.null(batch_id) | missing(batch_id)) {
    if ("WORKING_BATCH_ID" %in% ls(envir = neonmicrobe_env)) {
      # Save and reset batch environment
      saveParameters(verbose=FALSE)
      rm(list=ls(envir=batch_env), envir=batch_env)

      # Remove batch indicator variables
      prev_batch_id <- get("WORKING_BATCH_ID", envir = neonmicrobe_env)
      rm("WORKING_BATCH_ID", envir = neonmicrobe_env)
      rm("WORKING_BATCH_DIR", envir = neonmicrobe_env)

      if(verbose) message("Stepped outside of processing batch structure. No working batch ID is currently set. To undo this, run setBatch('", prev_batch_id, "')")
    }

    # If batch_id is provided, set WORKING_BATCH_ID to batch_id and load batch-specific parameters
  } else {
    # If already in a batch, save and reset batch environment (unless batch ID is the same)
    if ("WORKING_BATCH_ID" %in% ls(envir = neonmicrobe_env)) {
      if(identical(get("WORKING_BATCH_ID", envir=neonmicrobe_env), batch_id)) {
        if(verbose) message("Already in batch ID '", batch_id, "'.")
        return(invisible(NULL))
      }
      # Save and reset batch environment
      saveParameters(verbose=FALSE)
      rm(list=ls(envir=batch_env), envir=batch_env)
    }
    if(!dir.exists(batches_dir)) dir.create(batches_dir, recursive=TRUE)
    this_batch_dir <- file.path(batches_dir, batch_id)
    if(!dir.exists(this_batch_dir)) {
      message("Batch does not exist. Create a new batch using newBatch().")
      return(invisible(NULL))
    } else {
      assign("WORKING_BATCH_ID", batch_id, env=neonmicrobe_env)
      assign("WORKING_BATCH_DIR", normalizePath(this_batch_dir), env=neonmicrobe_env)

      # Load new batch_environment
      if(!identical(suppress_load_parameters, TRUE)) {
        loadParameters(file.path(this_batch_dir, "params.R"), batch_env, verbose=FALSE)
      }
      if(verbose) message("Switched to batch at ", this_batch_dir, ". Now working with processing batch ID '", batch_id, "'.")
    }
  }
}


#' Get Processing Batch Parameter(s)
#'
#' Gets parameter(s) associated with the currently loaded processing batch. Based on \code{\link[dada2]{getDadaOpt}}.
#'
#' @param params Parameters to get. By default (NULL), returns all parameters.
#'
#' @return Named list of parameters if multiple parameters are retrieved; else, single element (can be char, int, double).
#' @export
#'
#' @seealso \code{\link{setBatchParameter}}
#'
#' @examples
getBatchParameter <- function(params = NULL) {
  if ("WORKING_BATCH_ID" %in% ls(envir = neonmicrobe_env)) {
    if(is.null(params)) params <- ls(envir = batch_env)
    if(!all(params %in% ls(envir = batch_env))) {
      warning("Tried to get a nonexistent neonMicrobe processing batch parameter: ", params[!(params %in% ls(envir = batch_env))])
      params <- params[params %in% ls(envir = batch_env)]
    }
    params_select <- lapply(params, function(x) get(x, envir=batch_env))
    names(params_select) <- params
    if(length(params_select) == 1) params_select <- params_select[[1]]  # If just one param requested, return it alone
    return(params_select)
  } else {
    warning("Processing batch is not currently set. Use setBatch() to do so.")
    return(invisible(NULL))
  }
}


#' Set Processing Batch Parameter(s)
#'
#' Sets parameter(s) associated with the currently loaded processing batch. Based on \code{\link[dada2]{setDadaOpt}}.
#'
#' @param ... The processing parameters to set, along with their new value.
#' @param verbose Default TRUE. Whether to print a message confirming the assignment of a parameter.
#'
#' @return No value is returned.
#' @export
#'
#' @seealso \code{\link[dada2]{filterAndTrim}}, \code{\link[dada2]{setDadaOpt}}
#'
#' @details
#'
#' Quality filtering parameters, used as defaults in \code{\link[dada2]{filterAndTrim}}:
#'
#' MAX_EE_FWD, MAX_EE_REV: max. allowable expected errors in forward/reverse reads that pass filter
#' MIN_LEN_FWD, MIN_LEN_REV: min. allowable length of reads that pass filter in forward/reverse reads
#' TRUNC_Q_FWD, TRUNC_Q_REV: base quality score after which to truncate the sequence in forward/reverse reads
#' TRUNC_LEN_FWD, TRUNC_LEN_REV: length of read after which to truncate the sequence in forward/reverse reads. Reads shorter than this value are discarded.
#'
#' Parameters for DADA options, including heuristics: see \code{\link[dada2]{setDadaOpt}}
setBatchParameter <- function(..., verbose=TRUE) {
  if ("WORKING_BATCH_ID" %in% ls(envir = neonmicrobe_env)) {
    params <- getBatchParameter()
    args <- list(...)
    if(length(args)==1 && is.list(args[[1]])) {
      args <- args[[1]]
    }
    for(p_nm in names(args)) {
      # if(p_nm %in% names(params)) {
          # if(class(getBatchParameter(p_nm)) == class(args[[p_nm]])) {
            assign(p_nm, args[[p_nm]], envir = batch_env)
            if(verbose) message("Assigned to parameter '", p_nm, "': ", args[[p_nm]])

          # }
      # } else {
      #   warning(p_nm, " is not a valid processing batch parameter.")
      # }
    }
    # Save parameter changes to params file
    saveParameters(verbose=FALSE)
  } else {
    warning("Processing batch is not currently set. Use setBatch() to do so.")
    return(invisible(NULL))
  }
}


#' List Processing Batches
#'
#' Lists all processing batches available in the specified batches directory, along with basic information about them.
#'
#' @param batches_dir
#'
#' @return
#' @export
#'
#' @examples
listBatches <- function(batches_dir = file.path(NEONMICROBE_DIR_OUTPUTS, "batches")) {
  if(dir.exists(batches_dir)) {
    message("List of processing batches in ", batches_dir, ":")
    dirs <- list.dirs(batches_dir, full.names=TRUE, recursive=FALSE)
    about <- t(sapply(dirs, function(x) readLines(con = file.path(x, "ABOUT"), 2)))
    about_parsed <- sub(".*: ", "", about)
    about_parsed[,3] <- basename(about_parsed[,3])
    return(data.frame(
      ID = basename(dirs),
      created = about_parsed[,1],
      pkgVersion = about_parsed[,2],
      sequence_metadata = about_parsed[,3],
      row.names = NULL
    ))
  } else {
    message("Specified directory does not exist in specified location. Check 'batches_dir' argument, or create your first processing batch using newBatch().")
  }
}



#' Check Function Arguments against Batch-Specific Parameters
#'
#' If a processing batch is currently set, inserting this function into another function from
#' neonMicrobe will check the parent function's arguments against the parameters associated
#' with the current processing batch. If priority = "batch", the current batch's parameters would
#' take precedence.
#'
#' @details
#' If priority = "batch", then rather than using quality filtering parameters
#' defined on-the-fly as function arguments, the parent function would instead
#' use the quality filtering parameters associated with the current batch.
#' If priority = "batch", then rather than using the default output directory,
#' the parent function would instead use the batch-specific outputs.
#'
#' @param priority Must be "batch" or "arguments". Whether to defer to batch-specific processing parameters or function arguments when they conflict.
#' @param ... Format [function argument (string)] = [batch parameter (string)]. Arguments with a corresponding batch parameter, set to whatever the name of that batch parameter is, e.g. `"maxEE" = "MAX_EE_FWD`. The arguments will be checked against the batch parameters.
#' @param verbose Default TRUE. Whether to print a warning if the function argument and batch parameter do not match.
#' @param warn_no_batch Default FALSE. Whether to print a warning if there is no current processing batch.
#'
#' @return No value is returned
#'
#' @examples
checkArgsAgainstBatchParams <- function(priority = c("batch", "arguments"), ... , verbose = TRUE, warn_no_batch = FALSE) {
  priority <- match.arg(priority)
  dots <- list(...)
  if ("WORKING_BATCH_ID" %in% ls(envir = neonmicrobe_env)) {
    for(arg_nm in names(dots)) {
      par_val <- get(dots[[arg_nm]], envir=batch_env)
      arg_val <- get(arg_nm, envir=parent.frame())
      # If the batch parameter and the function argument DO NOT match...
      if(!identical(par_val, arg_val)) {
        if(identical(priority, "batch")) {
          par_val_print <- ifelse(is.character(par_val), shQuote(par_val), par_val)
          if(verbose) {
            warning("Argument '", arg_nm, "' does not match with batch parameter '", dots[[arg_nm]],
                    "'. Using batch parameter value: ", par_val_print, ". (Change behavior with 'priority' arg.)")
          }
          assign(arg_nm, par_val, envir=parent.frame())
        } else {
          arg_val_print <- ifelse(is.character(arg_val), shQuote(arg_val), arg_val)
          if(verbose) {
            warning("Argument '", arg_nm, "' does not match with batch parameter '", dots[[arg_nm]],
                    "'. Keeping function argument value: ", arg_val_print, ". (Change behavior with 'priority' arg.)")
          }
        }
      # If the function argument and the batch parameter DO match...
      } else {
        next
      }
    }
  } else {
    if(warn_no_batch) warning("Processing batch is not currently set. Use setBatch() to do so.")
  }
}

