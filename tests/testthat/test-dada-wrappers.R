
test_that("fastq in extdata are loadable", {
  fl_nm <- c("BMI_Plate23WellE7_16S_BDNB6_R1.fastq.gz",
             "BMI_Plate23WellE7_16S_BDNB6_R2.fastq.gz",
             "BMI_Plate37WellA12_16S_BJ8RK_R1.fastq.gz",
             "BMI_Plate37WellA12_16S_BJ8RK_R2.fastq.gz")
  files <- system.file("extdata", fl_nm, package="neonMicrobe")
  expect_equal(length(files), 4)
  expect_true(all(file.exists(files)))

  # put files in raw sequence directory for next tests
  setBaseDirectory(dirname(dirname(getwd())))
  file.copy(files, file.path(NEONMICROBE_DIR_SEQUENCE(), "16S", fl_nm), overwrite=FALSE)
})

test_that("trimPrimers16S does not fail", {
  data("seqmeta_greatplains")
  expect_error(trimPrimers16S(fn = fl_nm, in_subdir = "raw", out_subdir = "1_trimmed", meta = seqmeta_greatplains),
               NA)
})

test_that("qualityFilter16S does not fail", {
  expect_error(qualityFilter16S(fn = fl_nm, in_subdir = "1_trimmed", out_subdir = "2_filtered", meta = seqmeta_greatplains),
               NA)
})

# This test takes a long time (~15 min)
test_that("runDada16S does not fail", {
  expect_error(runDada16S(fn = fl_nm, in_subdir = "2_filtered", meta = seqmeta_greatplains, multithread = 8),
               NA)
})

test_that("no error if zero files are input", {
  # This may happen if, for instance, the user processes sequences
  # on a run-by-run basis without checking first that each run
  # has files within the set of input filenames.
  trimPrimers16S(fn = character(0), in_subdir = "raw", out_subdir = "1_trimmed", meta = seqmeta_greatplains)
  expect_error(trimPrimers16S(fn = character(0), in_subdir = "raw", out_subdir = "1_trimmed", meta = seqmeta_greatplains),
               NA)

})

# test_that("warning if meta does not cover all input files", {
#
# })


# test_fn_inner <- function(x) {
#   print(missing(x))
#   print(is.null(x))
#   print(x)
# }
# test_fn_outer <- function(y) test_fn_inner(y)
# test_fn_outer()
