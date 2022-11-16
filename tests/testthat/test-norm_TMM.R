library(rnaDualSeq)

setwd("./../..")

test_that("Right type of output", {
  host <- read_csv("data/Host.data.csv")
  phen <- read_csv("data/Pheno.csv")
  test_norm <- norm_TMM(host, phen)
  expect_that(is.matrix(test_norm), equals(TRUE) )
})
#
# test_that("Data files not valid", {
#   # test_read <- read_csv("fake/path")
#   expect_error(read_csv("fake/path"), "Not valid input!" )
# })
