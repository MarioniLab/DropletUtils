# This tests the cleanHashDrops() function.
# library(testthat); library(DropletUtils); source("test-cleanhash.R")

x <- rbind(
    rpois(1000, rep(c(100, 10), c(100, 900))),
    rpois(1000, rep(c(20, 100, 20), c(100, 100, 800))),
    rpois(1000, rep(c(30, 100, 30), c(200, 700, 100)))
)

set.seed(1000001)
test_that("cleanHashDrops works as expected", {
    x <- cbind(0, x, 1000)

    df <- cleanHashDrops(x)
    expect_identical(which(df$zero.ambient), 1L)
    expect_true(ncol(x) %in% which(df$high.ambient))
})

set.seed(1000001)
test_that("cleanHashDrops works correctly with controls", {
    x2 <- rbind(x, Y=c(0, rep(100, ncol(x) - 1)))
    df <- cleanHashDrops(x2, controls="Y")
    expect_identical(df$sum.controls, x2["Y",])
    expect_identical(which(df$low.controls), 1L)

    amb <- 1:4 
    names(amb) <- rownames(x2)
    df <- cleanHashDrops(x2, controls="Y", ambient=amb, contrib.method="control")
    expect_identical(df$ambient.scale[1], 0)
    expect_true(all(df$ambient.scale[-1]==100/4))
})
