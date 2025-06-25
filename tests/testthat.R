# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(strainspy)

library(R.utils) 
# We annoyingly need R.utils to be able to read a zipped csv file using read.csv
# Ignore note about import not being used, github workflow fails if this is removed

test_check("strainspy")
