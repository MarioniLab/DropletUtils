# This tests that the C++ code for emptyDrops does what it says it does.
# library(DropletUtils); library(testthat); source("test-emptydrops.R")

test_that("p-value calculations are correct", {
    # Simulating some counts.
    SIMSTUFF <- function(nbarcodes, ngenes) { 
        totals <- sample(1000, nbarcodes, replace=TRUE)
        ambient <- runif(ngenes)
        probs <- numeric(length(totals))
        for (i in seq_along(totals)) { 
            sim <- rmultinom(1, totals[i], prob=ambient)
            probs[i] <- sum(sim * log(ambient) - lfactorial(sim))
        }
        return(list(totals=totals, probs=probs, ambient=ambient))
    }

    # Setting up the reference function.
    REFFUN <- function(totals, probs, ambient, iter) {
        cumamb <- cumsum(ambient)
        logamb <- log(ambient)
        sumamb <- sum(ambient)

        n.above <- integer(length(totals))
        for (it in seq_len(iter)) {
            vec <- runif(max(totals), 0, sumamb)
            chosen <- findInterval(vec, cumamb)+1L
            
            for (i in seq_along(totals)) {
                curchosen <- chosen[seq_len(totals[i])]
                curprof <- tabulate(curchosen, nbins=length(ambient))
                curP <- sum(logamb * curprof - lfactorial(curprof))
                if (curP <= probs[i]) { 
                    n.above[i] <- n.above[i] + 1L
                }
            }
        }

        return(n.above)
    }

    # Checking the equivalence of our two implementations.
    set.seed(0)
    sim <- SIMSTUFF(100, 20)
    totals <- sim$totals
    probs <- sim$probs
    ambient.prof <- sim$ambient

    set.seed(100)
    ref <- REFFUN(totals, probs, ambient.prof, iter=100)
    set.seed(100)
    stats <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=100)
    expect_identical(stats, ref)

    # More barcodes (especially checking what happens with two barcodes of the same total).
    set.seed(0)
    sim <- SIMSTUFF(1001, 20)
    totals <- sim$totals
    probs <- sim$probs
    ambient.prof <- sim$ambient

    set.seed(100)
    ref <- REFFUN(totals, probs, ambient.prof, iter=100)
    set.seed(100)
    stats <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=100)
    expect_identical(stats, ref)

    # More genes.
    set.seed(0)
    sim <- SIMSTUFF(100, 200)
    totals <- sim$totals
    probs <- sim$probs
    ambient.prof <- sim$ambient

    set.seed(100)
    ref <- REFFUN(totals, probs, ambient.prof, iter=100)
    set.seed(100)
    stats <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=100)
    expect_identical(stats, ref)

    # More iterations.
    set.seed(0)
    sim <- SIMSTUFF(100, 20)
    totals <- sim$totals
    probs <- sim$probs
    ambient.prof <- sim$ambient

    set.seed(100)
    ref <- REFFUN(totals, probs, ambient.prof, iter=500)
    set.seed(100)
    stats <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=500)
    expect_identical(stats, ref)
})

test_that("emptyDrops runs to completion", {
    # Mocking up some counts.
    source(system.file("scripts", "mock_empty.R", package="DropletUtils"))
    limit <- 100
    e.out <- emptyDrops(my.counts, lower=limit)
    
    totals <- colSums(my.counts)
    expect_identical(as.integer(totals), e.out$Total)
    expect_true(all(is.na(e.out$LogProb[totals<=limit])))
    expect_true(all(!is.na(e.out$LogProb[totals>limit])))

    # Checking that the log-probability calculation is correct.
    ambient.prop <- edgeR::goodTuringProportions(Matrix::rowSums(my.counts[,totals<= limit]))
    valid <- which(totals > limit)
    collected <- numeric(nrow(e.out))
    for (v in valid) { 
        collected[v] <- dmultinom(my.counts[,v], size=totals[v], prob=ambient.prop, log=TRUE)
    }
    expect_equal(e.out$LogProb[valid], collected[valid])

    # Testing barcode ranks.    
    brout <- barcodeRanks(my.counts, lower=limit)
    expect_identical(brout$total, totals)
    expect_identical(brout$rank, rank(-totals, ties.method="average"))
    expect_true(all(is.na(brout$fitted[totals <= limit])))
    expect_true(all(!is.na(brout$fitted[totals > limit])))

    K <- brout$knee
    expect_true(all(e.out$FDR[totals >= K]==0))

    # Checking ambient tests.
    e.out2 <- emptyDrops(my.counts, test.args=list(test.ambient=TRUE))
    expect_identical(e.out$Total, e.out2$Total)
    expect_identical(e.out$Deviance[totals>limit], e.out2$Deviance[totals>limit])
    expect_true(all(!is.na(e.out2$LogProb[totals>0])))

    # Checking automatic retention options.
    e.out <- emptyDrops(my.counts, retain=K*0.6)
    expect_true(all(e.out$FDR[totals >= K*0.6]==0))
})
