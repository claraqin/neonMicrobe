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

test_that("seqmeta data can be loaded", {
  expect_silent(data("seqmeta_greatplains"))
})
