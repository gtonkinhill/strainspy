# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(strainspy)

if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils", repos = "https://cloud.r-project.org")
}

library(R.utils)

test_check("strainspy")
