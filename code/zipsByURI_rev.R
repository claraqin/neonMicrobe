# updated version of function, based on current version of function
# Lee F Stanish
# new arguments: data.format, dat, vars
# data.format: the format of the input data: can be the stackByTable format (same as before, no changes), an r data.frame containing the URLs,
# or a standalone .csv file of NEON data containing the URLs to download. 
# if data.format='rDataFrame'|'inputFile', then dat and vars are required. 
# dat: input data: either an r data.frame of NEON data with the target URLs to download, or a .csv file of NEON data with the targetURLs
# vars: name of the data field containing the URLs to download
# NOTE: pick.files values of TRUE will not work if data.format='rDataFrame' or 'inputFile'


zipsByUri <- function (filepath, 
                       data.format="stackByTable",
                       savepath = paste0(filepath, "/ECS_zipFiles"), 
                       pick.files = FALSE, 
                       check.size = TRUE, 
                       unzip = TRUE, 
                       saveZippedFiles = FALSE,
                       dat=NA,
                       vars=NA) {
  
  # Validate data.format #
  if(data.format == "rDataFrame") {
    if(nrow(dat)==0|is.na(vars)) {
      stop("dat and vars arguments required for specified data.format")
    } else {
      URLs.df <- dat
      URLsToDownload <- URLs.df[,which(names(URLs.df)==vars)]
      # verify URLs
      if(any(grepl("https", URLsToDownload)==FALSE)) {
        stop("at least one URL is not valid. Verify correct data and variable inputs")
      }
      if(length(URLsToDownload)==0) {
        stop("Unable to find valid URLs for input data")
      }
    }
  } else if(data.format == "inputFile") {
    URLs.df <- read.csv(dat, stringsAsFactors = FALSE)
    URLsToDownload <- URLs.df[,which(names(URLs.df)==vars)]
    # verify URLs
    if(any(grepl("https", URLsToDownload)==FALSE)) {
      stop("at least one URL is not valid. Verify correct data and variable inputs")
    }
    if(length(URLsToDownload)==0) {
      stop("Unable to find valid URLs for input data")
    }
  }
  
  
  #### Check for the variables file in the filepath
  if(data.format == "stackByTable") {
    files <- list.files(filepath, pattern = "variables")
    if (length(files) == 0) {
      stop("Variables file is not present in specified filepath.")
    }
    if (length(files) > 1) {
      stop("More than one variables file found in filepath.")
    }
    variablesFile <- utils::read.csv(paste(filepath, files, sep = "/"), 
                                     stringsAsFactors = FALSE)
    URLs <- variablesFile[variablesFile$dataType == "uri", ]
    allTables <- unique(URLs$table)
    URLsToDownload <- NA
    URLsNotToDownload <- NA
    if (pick.files == TRUE) {
      for (i in seq(along = allTables)) {
        suppressWarnings(tableData <- try(utils::read.csv(paste(filepath, 
                                                                "/", allTables[i], ".csv", sep = ""), stringsAsFactors = FALSE), 
                                          silent = TRUE))
        if (!is.null(attr(tableData, "class")) && attr(tableData, 
                                                       "class") == "try-error") {
          cat("Unable to find data for table:", allTables[i], 
              "\n")
          next
        }
        URLsPerTable <- names(tableData)[names(tableData) %in% 
                                           URLs$fieldName]
        for (j in URLsPerTable) {
          if (j %in% URLsToDownload | j %in% URLsNotToDownload) {
            next
          }
          resp <- readline(paste("Continuing will download", 
                                 length(tableData[, j]), "files for", j, "in", 
                                 allTables[i], "table. Do you want to include y/n: ", 
                                 sep = " "))
          if (resp %in% c("y", "Y")) {
            URLsToDownload <- c(URLsToDownload, tableData[, 
                                                          j])
          }
          else {
            URLsNotToDownload <- c(URLsNotToDownload, tableData[, 
                                                                j])
          }
        }
      }
    }
    else {
      for (i in seq(along = allTables)) {
        suppressWarnings(tableData <- try(utils::read.csv(paste(filepath, 
                                                                "/", allTables[i], ".csv", sep = ""), stringsAsFactors = FALSE), 
                                          silent = TRUE))
        if (!is.null(attr(tableData, "class")) && attr(tableData, 
                                                       "class") == "try-error") {
          cat("Unable to find data for table:", allTables[i], 
              "\n")
          next
        }
        URLsPerTable <- which(names(tableData) %in% URLs$fieldName)
        URLsToDownload <- c(URLsToDownload, unlist(tableData[, 
                                                             URLsPerTable]))
      }
      URLsToDownload <- unique(URLsToDownload)
    }
  }
  
  #Remove NA values from the list of URLs
  URLsToDownload <- URLsToDownload[!is.na(URLsToDownload)]
  if (length(URLsToDownload) == 0) {
    stop("There are no URLs other than NA for the stacked data.")
  }
  #Create the directory only if it doesn't already exist
  if (!dir.exists(savepath)) {
    dir.create(savepath)
  }
  #Loop to check existence and cumulative size of files
  cat("checking file sizes...\n")
  fileSize <- rep(NA, length(URLsToDownload))
  idx <- 0
  idxrem <- NA
  for (i in URLsToDownload) {
    idx <- idx + 1
    response <- httr::HEAD(i)
    if (is.null(httr::headers(response)[["Content-Length"]])) {
      cat(paste("No files found for url ", i, "\n", sep = ""))
      idxrem <- c(idxrem, idx)
    }
    else {
      fileSize[idx] <- as.numeric(httr::headers(response)[["Content-Length"]])
    }
  }
  totalFileSize <- convByteSize(sum(fileSize, na.rm = TRUE))
  if (check.size == TRUE) {
    resp <- readline(paste("Continuing will download", length(URLsToDownload), 
                           "files totaling approximately", totalFileSize, ". Do you want to proceed y/n: ", 
                           sep = " "))
    if (!(resp %in% c("y", "Y"))) 
      stop()
  }
  else {
    cat("Downloading", length(URLsToDownload), "files totaling approximately", 
        totalFileSize, ".\n")
  }
  # remove URLs with no data
  idxrem <- idxrem[-1]
  if (length(idxrem) > 0) {
    URLsToDownload <- URLsToDownload[-idxrem]
  }
  # copy zip files into folder
  numDownloads <- 0
  pb <- utils::txtProgressBar(style = 3)
  utils::setTxtProgressBar(pb, 1/(length(URLsToDownload) - 
                                    1))
  for (i in URLsToDownload) {
    dl <- try(downloader::download(i, paste(savepath, gsub("^.*\\/", 
                                                           "", i), sep = "/"), quiet = TRUE, mode = "wb"))
    if (!is.null(attr(dl, "class")) && attr(dl, "class") == 
        "try-error") {
      cat("Unable to download data for URL:", i, "\n")
      next
    }
    numDownloads <- numDownloads + 1
    utils::setTxtProgressBar(pb, numDownloads/(length(URLsToDownload) - 
                                                 1))
    if (unzip == TRUE && grepl("\\.zip|\\.ZIP", i)) {
      utils::unzip(paste(savepath, gsub("^.*\\/", "", i), 
                         sep = "/"), exdir = paste(savepath, gsub("^.*\\/|\\..*$", 
                                                                  "", i), sep = "/"), overwrite = TRUE)
      if (!saveZippedFiles) {
        unlink(paste(savepath, gsub("^.*\\/", "", i), 
                     sep = "/"), recursive = FALSE)
      }
    }
    else if (unzip == TRUE && grepl("\\.tar\\.gz", i)) {
      utils::untar(paste(savepath, gsub("^.*\\/", "", i), 
                         sep = "/"), exdir = paste(savepath, gsub("^.*\\/|\\..*$", 
                                                                  "", i), sep = "/"))
      if (!saveZippedFiles) {
        unlink(paste(savepath, gsub("^.*\\/", "", i), 
                     sep = "/"), recursive = FALSE)
      }
    }
    else if (unzip == TRUE && (grepl("\\.fastq\\.gz", i))) {
      if (!requireNamespace("R.utils", quietly = T)) {
        stop("Package R.utils is required for this function to work on fastq files. Please install and try again.")
      }
      R.utils::gunzip(paste(savepath, gsub("^.*\\/", "", 
                                           i), sep = "/"), remove = FALSE)
      if (!saveZippedFiles) {
        unlink(paste(savepath, gsub("^.*\\/", "", i), 
                     sep = "/"), recursive = FALSE)
      }
    }
    else if (grepl("\\.csv|\\.CSV", i)) {
      next
    }
    else if (unzip == TRUE && !(grepl("\\.zip|\\.ZIP", i) | 
                                grepl("\\.tar\\.gz", i))) {
      cat("Unable to unzip data for URL:", i, "\n")
    }
  }
  utils::setTxtProgressBar(pb, 1)
  close(pb)
  cat(numDownloads, "file(s) successfully downloaded to", savepath, 
      "\n", sep = " ")
}