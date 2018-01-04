# Test that function runs to completion.
# library(DropletUtils); library(testthat); source("test-default.R")

test_that("defaultDrops runs to completion", {
    # Mock up counts
    source(system.file("scripts", "mock_empty.R", package="DropletUtils"))
    out <- defaultDrops(my.counts)
   
    # Should always call at least one cell (100th %ile cell)
    expect_true(sum(out)>0)
})
