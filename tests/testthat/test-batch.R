test_that("warning if base dir not set", {
  expect_warning(NEONMICROBE_DIR_BASE())
  setBaseDirectory()
})

test_that("batch_env exists", {
  expect_true(exists("batch_env"))
})

test_that("batch_env persists", {
  newBatch(system.file("extdata", "seqmeta_greatplains.csv", package="neonMicrobe"), batch_id="test", overwrite = TRUE)
  expect_true(exists("batch_env"))
  setBatch()
  expect_true(exists("batch_env"))

  #cleanup
  setBatch()
  unlink(file.path(NEONMICROBE_DIR_BATCHES(), "test"), recursive=TRUE)
})

test_that("batch indicator variables work", {
  setBatch()
  expect_false("WORKING_BATCH_ID" %in% ls(envir = neonmicrobe_env))
  expect_false("WORKING_BATCH_DIR" %in% ls(envir = neonmicrobe_env))
  newBatch(system.file("extdata", "seqmeta_greatplains.csv", package="neonMicrobe"), batch_id="test", overwrite = TRUE)
  expect_true("WORKING_BATCH_ID" %in% ls(envir = neonmicrobe_env))
  expect_true("WORKING_BATCH_DIR" %in% ls(envir = neonmicrobe_env))

  #cleanup
  setBatch()
  unlink(file.path(NEONMICROBE_DIR_BATCHES(), "test"), recursive=TRUE)
})

test_that("params cannot be set outside a batch", {
  setBatch()
  expect_warning(setBatchParam(TEST="test"))
})

test_that("params can be set within a batch", {
  newBatch(system.file("extdata", "seqmeta_greatplains.csv", package="neonMicrobe"), batch_id="test", overwrite = TRUE)
  setBatchParam(TEST="test")
  expect_true("TEST" %in% names(getBatchParam()))

  #cleanup
  setBatch()
  unlink(file.path(NEONMICROBE_DIR_BATCHES(), "test"), recursive=TRUE)
})

test_that("setting param within a batch does not affect params outside batch", {
  newBatch(system.file("extdata", "seqmeta_greatplains.csv", package="neonMicrobe"), batch_id="test", overwrite = TRUE)
  setBatchParam(TEST="test")
  expect_false("TEST" %in% ls(neonmicrobe_env))
  expect_true(exists("TEST", batch_env))
  setBatch()
  expect_false(exists("TEST", batch_env))
  expect_equal(ls(batch_env), character(0))

  #cleanup
  setBatch()
  unlink(file.path(NEONMICROBE_DIR_BATCHES(), "test"), recursive=TRUE)
})

test_that("batch params can be saved to file", {
  setBatch() # Ensures that we start outside of a batch
  expect_warning(saveParams(verbose=FALSE))
  batch_id1 <- paste0("test_", rnorm(1))
  batch_id2 <- paste0("test_", rnorm(1))
  newBatch(system.file("extdata", "seqmeta_greatplains.csv", package="neonMicrobe"), batch_id=batch_id1, overwrite = TRUE)
  setBatchParam(TEST="test")
  newBatch(system.file("extdata", "seqmeta_greatplains.csv", package="neonMicrobe"), batch_id=batch_id2, overwrite = TRUE)
  expect_false(exists("TEST", batch_env))
  loadParams(file.path(NEONMICROBE_DIR_BATCHES(), batch_id1, "params.R"), verbose=FALSE)
  expect_true(exists("TEST", batch_env))

  #cleanup
  setBatch()
  unlink(file.path(NEONMICROBE_DIR_BATCHES(), batch_id1), recursive=TRUE)
  unlink(file.path(NEONMICROBE_DIR_BATCHES(), batch_id2), recursive=TRUE)
})

test_that("checkArgsAgainstBatchParams can change function args", {
  newBatch(system.file("extdata", "seqmeta_greatplains.csv", package="neonMicrobe"), batch_id="test", overwrite = TRUE)
  setBatchParam(TEST = "parameter_value")
  test_fn <- function(x, prioritize_params=TRUE) {
    original_val <- x
    if(prioritize_params) {
      checkArgsAgainstBatchParams("x" = "TEST", priority="batch")
    } else {
      checkArgsAgainstBatchParams("x" = "TEST", priority="arguments")
    }
    new_val <- x
    return(c(original_val, new_val))
  }
  expect_warning(out <- test_fn("original_value"))
  expect_identical(out[1], "original_value")
  expect_identical(out[2], "parameter_value")

  #cleanup
  setBatch()
  unlink(file.path(NEONMICROBE_DIR_BATCHES(), "test"), recursive=TRUE)
})

test_that("checkArgsAgainstBatchParams can leave function args intact if desired", {
  newBatch(system.file("extdata", "seqmeta_greatplains.csv", package="neonMicrobe"), batch_id="test", overwrite = TRUE)
  setBatchParam(TEST = "parameter_value")
  test_fn <- function(x, prioritize_params=TRUE) {
    original_val <- x
    if(prioritize_params) {
      checkArgsAgainstBatchParams("x" = "TEST", priority="batch")
    } else {
      checkArgsAgainstBatchParams("x" = "TEST", priority="arguments")
    }
    new_val <- x
    return(c(original_val, new_val))
  }
  expect_warning(out_warnonly <- test_fn("original_value", prioritize_params=FALSE))
  expect_identical(out_warnonly[1], out_warnonly[2])
  expect_identical(out_warnonly[1], "original_value")

  #cleanup
  setBatch()
  unlink(file.path(NEONMICROBE_DIR_BATCHES(), "test"), recursive=TRUE)
})

test_that("checkArgsAgainstBatchParams can change variables defined WITHIN function", {
  newBatch(system.file("extdata", "seqmeta_greatplains.csv", package="neonMicrobe"), batch_id="test", overwrite = TRUE)
  setBatchParam(TEST = "parameter_value")
  test_fn <- function(x) {
    original_val <- changeable_val <- x
    checkArgsAgainstBatchParams("changeable_val" = "TEST")
    new_val <- changeable_val
    return(c(original_val, new_val))
  }
  expect_warning(out <- test_fn("original_val"))
  expect_false(identical(out[1], out[2]))

  #cleanup
  setBatch()
  unlink(file.path(NEONMICROBE_DIR_BATCHES(), "test"), recursive=TRUE)
})

test_that("left side of checkArgsAgainstBatchParams args do not have to be enclosed in quotes", {
  newBatch(system.file("extdata", "seqmeta_greatplains.csv", package="neonMicrobe"), batch_id="test", overwrite = TRUE)
  setBatchParam(TEST = "parameter_value")
  test_fn <- function(x) {
    original_val <- x
    checkArgsAgainstBatchParams(x = "TEST")
    new_val <- x
    return(c(original_val, new_val))
  }
  expect_warning(out <- test_fn("original_val"))
  expect_false(identical(out[1], out[2]))
  expect_equal(out[2], "parameter_value")

  #cleanup
  setBatch()
  unlink(file.path(NEONMICROBE_DIR_BATCHES(), "test"), recursive=TRUE)
})

test_that("checkArgsAgainstBatchParams is tolerant to nonexistent arg", {

})

test_that("outputs directory changes with batch", {
  setBatch()
  dir1 <- NEONMICROBE_DIR_OUTPUTS()
  newBatch(system.file("extdata", "seqmeta_greatplains.csv", package="neonMicrobe"), batch_id="test", overwrite=TRUE)
  dir2 <- NEONMICROBE_DIR_OUTPUTS()
  expect_false(identical(dir1, dir2))

  #cleanup
  setBatch()
  unlink(file.path(NEONMICROBE_DIR_BATCHES(), "test"), recursive=TRUE)
})

test_that("first new batch does not throw error", {
  setBaseDirectory()
  batch_id <- paste0("test_", rnorm(1))
  expect_error(newBatch(system.file("extdata", "seqmeta_greatplains.csv", package="neonMicrobe"), batch_id=batch_id, overwrite=TRUE), NA)

  #cleanup
  setBatch()
  unlink(file.path(NEONMICROBE_DIR_BATCHES(), batch_id), recursive=TRUE)
})

test_that("Final cleanup!", {
  expect_equal(length(dir(NEONMICROBE_DIR_BATCHES())), 0)
  if(length(dir(NEONMICROBE_DIR_BATCHES())) == 0) {
    unlink(NEONMICROBE_DIR_BATCHES(), recursive=TRUE)
  }
})
