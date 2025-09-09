library('flowCore')

read_fcs <- function(filename) {
  if (isFCSfile(filename)) {
    read.FCS(filename)
  }
  else {
    stop("File is not of .fcs type")
  }
}