library(rnaDualSeq)

setwd("./../..")

test_that("Data files valid", {
  host <- read_csv("data/Host.data.csv")
  phen <- read_csv("data/Pheno.csv")
  test_norm <- norm_TMM(host, phen)
  test_DE <- identifyDE(test_norm, phen)
  expect_that(is.matrix(test_norm), equals(TRUE) )
})
