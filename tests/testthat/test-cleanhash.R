# This tests the cleanHashDrops() function.
# library(testthat); library(DropletUtils); source("test-cleanhash.R")

set.seed(1000001)
test_that("cleanHashDrops works as expected", {
    x <- rbind(
        rpois(1000, rep(c(100, 10), c(100, 900))),
        rpois(1000, rep(c(20, 100, 20), c(100, 100, 800))),
        rpois(1000, rep(c(30, 100, 30), c(200, 700, 100)))
    )
    x <- cbind(0, x, 1000)

    df <- cleanHashDrops(x)
    expect_identical(which(df$zero.ambient), 1L)
    expect_identical(which(df$high.ambient), ncol(x))

    # Works when we throw in some control features.
    x2 <- rbind(x, Y=c(0, rep(100, ncol(x) - 1)))
    df <- cleanHashDrops(x2, controls="Y")
    expect_identical(df$sum.controls, x2["Y",])
    expect_identical(which(df$low.controls), 1L)
})
