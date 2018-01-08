# This checks out the barcodeRanks and defaultDrops functions.
# library(DropletUtils); library(testthat); source("test-misc.R")

# Mocking up some counts.
set.seed(100)
my.counts <- DropletUtils:::simCounts()
totals <- Matrix::colSums(my.counts)

test_that("barcodeRanks runs to completion", {
    limit <- 100
    brout <- barcodeRanks(my.counts, lower=limit)
    expect_identical(brout$total, totals)
    expect_identical(brout$rank, rank(-totals, ties.method="average"))
    expect_true(all(is.na(brout$fitted[totals <= limit])))

    # Trying again with a higher limit.
    limit2 <- 200
    brout2 <- barcodeRanks(my.counts, lower=limit2)
    expect_identical(brout, brout2)

    # Specifying the boundaries.
    bounds <- c(200, 1000)
    brout3 <- barcodeRanks(my.counts, lower=limit, fit.bounds=bounds)
    is.okay <- totals > bounds[1] & totals < bounds[2]
    expect_true(all(is.na(brout3$fitted[!is.okay])))
    expect_true(all(!is.na(brout3$fitted[is.okay])))

    # Trying out silly inputs.
    expect_error(barcodeRanks(my.counts[,0]), "insufficient")
    expect_error(barcodeRanks(my.counts[0,]), "insufficient")
})

test_that("defaultDrops runs to completion", {
    out <- defaultDrops(my.counts)
   
    # Should always call at least one cell (100th %ile cell)
    expect_true(sum(out)>0)

    out <- defaultDrops(my.counts, lower.prop=0) # should keep all non-zero cells.
    expect_true(all(out | totals==0))

    out <- defaultDrops(my.counts, upper.quant=1, lower.prop=1) # as it's >, not >=.
    expect_true(!any(out))

    # Works alright on silly inputs.
    expect_identical(logical(0), defaultDrops(my.counts[,0]))
    expect_identical(logical(ncol(my.counts)), defaultDrops(my.counts[0,]))
})
