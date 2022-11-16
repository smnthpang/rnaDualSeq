library(rnaDualSeq)

setwd("./../..")

test_that("Data files valid", {
  test_read <- read_csv("data/Host.data.csv")
  expect_that(typeof(test_read), equals("list") )
})

test_that("Data files not valid", {
  # test_read <- read_csv("fake/path")
  expect_error(read_csv("fake/path"), "Not valid input!" )
})

