# This tests the maximum present ambience function.
# library(testthat); library(DropletUtils); source("setup.R"); source("test-amb-sparse.R")

# Making up some data.
set.seed(1100000)
ambient <- c(runif(900, 0, 0.1), runif(100))
y <- rpois(1000, ambient * 50)
y2 <- y + rpois(1000, 5) # adding some actual biology.

test_that("ambientContribSparse works correctly", {
    out <- ambientContribSparse(y, ambient=ambient, prop=0.5)
    expect_identical(unname(out), median(y/ambient))

    out <- ambientContribSparse(cbind(y, y2), ambient=ambient, prop=0.5)
    expect_identical(unname(out[1]), median(y/ambient))
    expect_identical(unname(out[2]), median(y2/ambient))

    # Handles zeroes in the ambient vector.
    ambient <- c(numeric(900), runif(100))
    out <- ambientContribSparse(cbind(y, y2), ambient=ambient, prop=0.5)
    keep <- 901:1000
    expect_identical(out, ambientContribSparse(cbind(y, y2)[keep,], ambient=ambient[keep], prop=0.5))
})
