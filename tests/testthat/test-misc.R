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

    # Respecting column names.
    alt <- my.counts
    colnames(alt) <- sprintf("BARCODE_%i", seq_len(ncol(alt)))
    brout2 <- barcodeRanks(alt)
    expect_identical(rownames(brout2), colnames(alt))
    expect_identical(names(brout2$rank), NULL)
    expect_identical(names(brout2$total), NULL)
    expect_identical(names(brout2$fitted), NULL)

    # Trying out silly inputs.
    expect_error(barcodeRanks(my.counts[,0]), "insufficient")
    expect_error(barcodeRanks(my.counts[0,]), "insufficient")
})

test_that("barcodeRanks' excluder works correctly", {
    brout <- barcodeRanks(my.counts)
    keep <- brout$total >= 100 & !duplicated(brout$total)
    x <- log10(brout$rank[keep])
    y <- log10(brout$total[keep])

    o <- order(x)
    x <- x[o]
    y <- y[o]

    # Compares correctly to a reference.
    edge.out <- DropletUtils:::.find_curve_bounds(x=x, y=y, exclude.from=100) 
    ref.out <- DropletUtils:::.find_curve_bounds(x=tail(x, -100), y=tail(y, -100), exclude.from=0) 
    expect_identical(edge.out, ref.out+100)

    edge.outx <- DropletUtils:::.find_curve_bounds(x=x, y=y, exclude.from=200) 
    ref.outx <- DropletUtils:::.find_curve_bounds(x=tail(x, -200), y=tail(y, -200), exclude.from=0) 
    expect_false(identical(edge.outx, ref.outx+200))

    # Proper edge behavior.
    edge.out2 <- DropletUtils:::.find_curve_bounds(x=x, y=y, exclude.from=0) 
    expect_identical(edge.out[2], edge.out2[2])
    expect_false(identical(edge.out[1], edge.out2[1]))

    edge.out3 <- DropletUtils:::.find_curve_bounds(x=x, y=y, exclude.from=Inf)
    expect_identical(unname(edge.out3[1]), length(y)-1)
    expect_identical(unname(edge.out3[2]), length(y)-1)

    # Works properly when put together. 
    ref <- barcodeRanks(my.counts)
    brout <- barcodeRanks(my.counts, exclude.from=100)
    expect_identical(ref, brout)

    brout2 <- barcodeRanks(my.counts, exclude.from=200)
    expect_false(identical(ref, brout2))

    brout3 <- barcodeRanks(my.counts, exclude.from=Inf)
    expect_false(identical(ref, brout2))
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
