# This tests the inferAmbience() function.
# library(testthat); library(DropletUtils); source("test-inferamb.R")

test_that("inferAmbience works correctly in the expected case", {
    mat <- rbind(c(1:100, 1000:1200))
    out <- inferAmbience(mat)
    expect_equal(out, 50.5)
})

test_that("inferAmbience works correctly with all-zero counts", {
    mat <- matrix(0, 10, 1000)
    rownames(mat) <- LETTERS[1:10]
    out <- inferAmbience(mat)
    expect_identical(out, mat[,1])
})
