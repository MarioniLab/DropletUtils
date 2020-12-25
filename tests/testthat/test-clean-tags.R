# This tests the cleanTagCounts() function.
# library(testthat); library(DropletUtils); source("test-clean-tags.R")

x <- rbind(
    rpois(1000, rep(c(100, 10), c(100, 900))),
    rpois(1000, rep(c(20, 100, 20), c(100, 100, 800))),
    rpois(1000, rep(c(30, 100, 30), c(200, 700, 100)))
)

test_that("cleanTagCounts works without controls", {
    x <- cbind(0, x, 1000)
    df <- cleanTagCounts(x)
    expect_identical(which(df$zero.ambient), 1L)
    expect_true(ncol(x) %in% which(df$high.ambient))
})

test_that("cleanTagCounts works correctly with controls", {
    x2 <- rbind(x, Y=c(1000, rep(100, ncol(x) - 1)))
    df <- cleanTagCounts(x2, controls="Y")
    expect_identical(df$sum.controls, x2["Y",])
    expect_identical(which(df$high.controls), 1L)
})

test_that("cleanTagCounts works correctly with exclusives", {
    x2 <- rbind(x, 
        Y=c(1000, 10, 1000, rep(10, ncol(x) - 3)),
        Z=c(1000, 1000, 10, rep(10, ncol(x) - 3))
    )
    df <- cleanTagCounts(x2, ambient=c(10, 20, 30, Y=10, Z=10), exclusive=c("Y", "Z"))
    expect_identical(df$ambient.scale[1:3], c(100, 1, 1))
    expect_true(1L %in% which(df$high.ambient))
})
