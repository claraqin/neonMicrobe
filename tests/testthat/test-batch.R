
# test_that("multiplication works", {
#   expect_equal(2 * 2, 4)
#   expect_equal(2 * 8, 16)
# })

test_that("warning if base dir not set", {
  expect_warning(NEONMICROBE_DIR_BASE())
})

test_that("getwd() in devtools::test() ends in tests/testthat/", {
  expect_identical(basename(getwd()), "testthat")
  expect_identical(basename(dirname(getwd())), "tests")
})

test_that("makeDataDirectories creates data directories in base dir", {
  setBaseDirectory(dirname(dirname(getwd())))
  makeDataDirectories(check_location=FALSE)
  expect_true(dir.exists(file.path(NEONMICROBE_DIR_BASE(), "data")))
  expect_true(dir.exists(file.path(NEONMICROBE_DIR_BASE(), "outputs")))
})

test_that("batch_env exists", {
  expect_true(exists("batch_env"))
})

test_that("batch_env persists", {
  newBatch(seq_meta_file=file.path(NEONMICROBE_DIR_BASE(), "/data/sequence_metadata/mmg_soilMetadata_ITS_2021-03-11080914.csv"), batch_id="z", set_batch = FALSE, overwrite = TRUE)
  expect_true(exists("batch_env"))
  setBatch()
  expect_true(exists("batch_env"))
})

test_that("batch indicator variables work", {
  setBatch()
  expect_false("WORKING_BATCH_ID" %in% ls(envir = neonmicrobe_env))
  expect_false("WORKING_BATCH_DIR" %in% ls(envir = neonmicrobe_env))
  newBatch(seq_meta_file=file.path(NEONMICROBE_DIR_BASE(), "/data/sequence_metadata/mmg_soilMetadata_ITS_2021-03-11080914.csv"), batch_id="z", set_batch = TRUE, overwrite = TRUE)
  expect_true("WORKING_BATCH_ID" %in% ls(envir = neonmicrobe_env))
  expect_true("WORKING_BATCH_DIR" %in% ls(envir = neonmicrobe_env))
})

test_that("params cannot be set outside a batch", {
  setBatch()
  expect_warning(setBatchParameter(TEST="test"))
})

test_that("params can be set within a batch", {
  newBatch(seq_meta_file=file.path(NEONMICROBE_DIR_BASE(), "/data/sequence_metadata/mmg_soilMetadata_ITS_2021-03-11080914.csv"), batch_id="z", set_batch = TRUE, overwrite = TRUE)
  setBatchParameter(TEST="test")
  expect_true("TEST" %in% names(getBatchParameter()))
})

test_that("setting param within a batch does not affect params outside batch", {
  newBatch(seq_meta_file=file.path(NEONMICROBE_DIR_BASE(), "/data/sequence_metadata/mmg_soilMetadata_ITS_2021-03-11080914.csv"), batch_id="z", set_batch = TRUE, overwrite = TRUE)
  setBatchParameter(TEST="test")
  expect_false("TEST" %in% ls(neonmicrobe_env))
  expect_true(exists("TEST", batch_env))
  setBatch()
  expect_false(exists("TEST", batch_env))
  expect_equal(ls(batch_env), character(0))
})

test_that("batch params can be saved to file", {
  setBatch() # Ensures that we start outside of a batch
  expect_warning(saveParameters(verbose=FALSE))
  newBatch(seq_meta_file=file.path(NEONMICROBE_DIR_BASE(), "/data/sequence_metadata/mmg_soilMetadata_ITS_2021-03-11080914.csv"), batch_id="z", set_batch = TRUE, overwrite = TRUE)
  setBatchParameter(TEST="test")
  expect_silent(saveParameters(verbose=FALSE))
  newBatch(seq_meta_file=file.path(NEONMICROBE_DIR_BASE(), "/data/sequence_metadata/mmg_soilMetadata_ITS_2021-03-11080914.csv"), batch_id="y", set_batch = TRUE, overwrite = TRUE)
  expect_false(exists("TEST", batch_env))
  loadParameters(file.path(NEONMICROBE_DIR_OUTPUTS(), "batches", "z", "params.R"), verbose=FALSE)
  expect_true(exists("TEST", batch_env))
})

test_that("checkArgsAgainstBatchParams works", {
  newBatch(seq_meta_file=file.path(NEONMICROBE_DIR_BASE(), "/data/sequence_metadata/mmg_soilMetadata_ITS_2021-03-11080914.csv"), batch_id="test", set_batch = TRUE, overwrite = TRUE)
  setBatchParameter(TEST = "parameter_value")
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

  expect_warning(out_warnonly <- test_fn("original_value", prioritize_params=FALSE))
  expect_identical(out_warnonly[2], "original_value")
})
