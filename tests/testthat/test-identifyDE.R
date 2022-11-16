library(rnaDualSeq)

setwd("./../..")

test_that("column names are valid", {
  host <- read_csv("data/Host.data.csv")
  phen <- read_csv("data/Pheno.csv")
  test_norm <- norm_TMM(host, phen)
  test_DE <- identifyDE(test_norm, phen)
  #colnames(test_DE[[1]])[1]
  expect_that(colnames(test_DE[[1]])[1], equals("logFC") )
  expect_that(colnames(test_DE[[1]])[2], equals("AveExpr") )
  expect_that(colnames(test_DE[[1]])[3], equals("t"  ) )
  expect_that(colnames(test_DE[[1]])[4], equals( "P.Value" ) )
  expect_that(colnames(test_DE[[1]])[5], equals("adj.P.Val") )
  expect_that(colnames(test_DE[[1]])[6], equals("B"      ) )
})
