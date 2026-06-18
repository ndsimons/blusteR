library(testthat)
# We source the package R files directly for testing without install
pkg_dir <- file.path(dirname(dirname(getwd())), "R")
if (dir.exists(pkg_dir)) {
  for (f in sort(list.files(pkg_dir, pattern = "\\.R$", full.names = TRUE))) {
    source(f)
  }
}
