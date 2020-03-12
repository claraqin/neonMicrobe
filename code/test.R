# test read fastq.gz directly from https
library(ShortRead)

# example code
fastqDir <- system.file("extdata", "E-MTAB-1147", package = "ShortRead")
fastqPath <- list.files(fastqDir, pattern = ".fastq.gz$", full = TRUE)[1]
reads <- readFastq(fastqPath)
reads

# download from AWS to temp and read locally
test.file <- "https://neon-microbial-raw-seq-files.s3.data.neonscience.org/2017/150925_M02149_0250_000000000-AHCGJ/BART.002.39.2.27.M.20140821.fastq.gz"
temp <- tempfile()
download.file(test.file, temp)
test.reads <- readFastq(temp)
sread(test.reads)
unlink(temp)

# directly read from https
# works for csv
data <- read.csv("http://apps.fs.fed.us/fiadb-downloads/CSV/LICHEN_SPECIES_SUMMARY.csv")
# doesn't work for fastq
readFastq(test.file)
readFastq(url(test.file))
readFastq(gzcon(url(test.file)))
readFastq(gzfile(url(test.file)))
