# This tests that the C++ code for emptyDrops does what it says it does.
# library(DropletUtils); library(testthat); source("test-empty.R")

test_that("multinomial calculations are correct", {
    set.seed(1000)
    prop <- runif(100)
    X <- matrix(rpois(100000, lambda=prop*2), nrow=length(prop))

    # Assumes that input has sum of 1!
    prop <- prop/sum(prop)

    # Checking it's the same regardless of the input matrix.
    dense.out <- DropletUtils:::.compute_multinom_prob(X, prop)
    sparse.out1 <- DropletUtils:::.compute_multinom_prob(as(X, "dgCMatrix"), prop)
    expect_equal(dense.out, sparse.out1)
    sparse.out2 <- DropletUtils:::.compute_multinom_prob(as(X, "dgTMatrix"), prop)
    expect_equal(dense.out, sparse.out2)

    # Checking it's the same compared to a reference.
    ref <- apply(X, 2, dmultinom, prob=prop, log=TRUE)
    expect_equal(ref, dense.out+lfactorial(colSums(X)))
})

test_that("p-value calculations are correct", {
    # Simulating some multinomial probabilities. 
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

    # Checking for uniformity:
#    set.seed(0)
#    sim <- SIMSTUFF(10000, 1000)
#    totals <- sim$totals
#    probs <- sim$probs
#    ambient.prof <- sim$ambient
#    set.seed(100)
#    stats <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=1000)
})

test_that("emptyDrops runs to completion", {
    # Mocking up some counts.
    set.seed(1000)
    my.counts <- DropletUtils:::simCounts()
    limit <- 100
    e.out <- emptyDrops(my.counts, lower=limit)
   
    totals <- Matrix::colSums(my.counts)
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

    # Checking ambient tests.
    e.out2 <- emptyDrops(my.counts, lower=limit, test.ambient=TRUE)
    expect_identical(e.out$Total, e.out2$Total)
    expect_identical(e.out$LogProb[totals>limit], e.out2$LogProb[totals>limit])
    expect_true(all(!is.na(e.out2$LogProb[totals>0])))

    # Checking the ignore argument works correctly.
    set.seed(1001)
    e.out3a <- emptyDrops(my.counts, lower=limit, ignore=200)
    e.out3a$FDR <- NULL
    set.seed(1001)
    survivors <- totals <= 100 | totals > 200
    new.counts <- my.counts[,survivors]
    e.out3b <- emptyDrops(new.counts, lower=limit, ignore=NULL)
    e.out3b$FDR <- NULL
    expect_equal(e.out3a[survivors,], e.out3b) 

    # Checking retention options.
    K <- barcodeRanks(my.counts, lower=limit)$knee
    expect_true(all(e.out$FDR[totals >= K]==0))
    expect_true(!all(e.out$FDR[totals < K]==0))

    e.outK <- emptyDrops(my.counts, retain=K*0.6)
    expect_true(all(e.outK$FDR[totals >= K*0.6]==0))
    expect_true(!all(e.outK$FDR[totals < K*0.6]==0))
})

