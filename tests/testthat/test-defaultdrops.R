#Test that function runs to completion
# library(DropletUtils); library(testthat); source("test-defaultdrops.R")

test_that("defaultDrops runs to completion", {
  #mock up counts
  source(system.file("scripts", "mock_empty.R", package="DropletUtils"))
  
  out = defaultDrops(my.counts)
  #should always call at least one cell (100th %ile cell)
  expect_true(sum(out)>0)
})
