#Test that function runs to completion
# library(DropletUtils); library(testthat); source("test-emptydrops.R")

test_that("cellrangerCall runs to completion", {
  #mock up counts
  source(system.file("scripts", "mock_empty.R", package="DropletUtils"))
  
  out = cellrangerCall(my.counts)
  #should always call at least one cell (100th %ile cell)
  expect_true(sum(out)>0)
})
