library(rnaDualSeq)

setwd("./../..")

test_that("New labels", {
  host <- read_csv("data/Host.data.csv")
  phen <- read_csv("data/Pheno.csv")
  test_norm <- norm_TMM(host, phen)
  test_DE <- identifyDE(test_norm, phen)
  volcanoPlot("host", test_DE, timeperiods = c("3h", "5h", "7h", "19h", "26h") ) # changed time periods for labelling

  expect_that(file.exists("Results"), equals(TRUE) )# Results folder is created
  expect_that(file.exists("Results/VolcanoPlots/host_3h.pdf"), equals(TRUE) )

})
