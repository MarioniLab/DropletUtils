# This checks out the barcodeRanks and defaultDrops functions.
# library(DropletUtils); library(testthat); source("test-misc.R")

# Mocking up some counts.
set.seed(100)
my.counts <- DropletUtils:::simCounts()
totals <- Matrix::colSums(my.counts)

test_that("barcodeRanks runs to completion", {
    limit <- 100
    brout <- barcodeRanks(my.counts, lower=limit)
    expect_equal(brout$total, totals)
    expect_identical(brout$rank, rank(-totals, ties.method="average"))

    # Trying again with a higher limit.
    limit2 <- 200
    brout2 <- barcodeRanks(my.counts, lower=limit2)
    expect_identical(brout, brout2)

    # Specifying the boundaries.
    bounds <- c(200, 1000)
    brout3 <- barcodeRanks(my.counts, lower=limit, fit.bounds=bounds)
    knee <- metadata(brout3)$knee
    expect_true(knee >= bounds[1] && knee <= bounds[2])

    # Respecting column names.
    alt <- my.counts
    colnames(alt) <- sprintf("BARCODE_%i", seq_len(ncol(alt)))
    brout2 <- barcodeRanks(alt)
    expect_identical(rownames(brout2), colnames(alt))
    expect_identical(names(brout2$rank), NULL)
    expect_identical(names(brout2$total), NULL)

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
