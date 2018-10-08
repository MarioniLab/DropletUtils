# This tests that the C++ code for emptyDrops does what it says it does.
# library(DropletUtils); library(testthat); source("test-empty.R")

test_that("multinomial calculations are correct without alpha", {
    set.seed(1000)
    prop <- runif(100)
    X <- matrix(rpois(100000, lambda=prop*2), nrow=length(prop))

    # DropletUtils assumes that proportions sum to unity!
    prop <- prop/sum(prop)

    # Checking it's the same regardless of the input matrix.
    dense.out <- DropletUtils:::.compute_multinom_prob_data(X, prop)
    sparse.out1 <- DropletUtils:::.compute_multinom_prob_data(as(X, "dgCMatrix"), prop)
    expect_equal(dense.out, sparse.out1)
    sparse.out2 <- DropletUtils:::.compute_multinom_prob_data(as(X, "dgTMatrix"), prop)
    expect_equal(dense.out, sparse.out2)

    # Checking it's the same compared to a reference.
    ref <- apply(X, 2, dmultinom, prob=prop, log=TRUE)
    combined <- dense.out + DropletUtils:::.compute_multinom_prob_rest(colSums(X))
    expect_equal(ref, combined)
})

ddirmultinom <- function(x, prob, alpha) 
# A helper function for computing the Dirichlet-multionial PDF.
{
    size <- sum(x)
    lfactorial(size) + lgamma(alpha) - lgamma(size + alpha) +
        sum(lgamma(x + prob * alpha) - lfactorial(x) - lgamma(prob * alpha))
} 

test_that("multinomial calculations are correct with alpha", {
    set.seed(2000)
    prop <- runif(100)
    X <- matrix(rpois(100000, lambda=prop*2), nrow=length(prop))
    prop <- prop/sum(prop)

    # Checking it's the same regardless of the input matrix.
    for (alpha in c(1, 10, 100)) { 
        dense.out <- DropletUtils:::.compute_multinom_prob_data(X, prop, alpha=alpha)
        sparse.out1 <- DropletUtils:::.compute_multinom_prob_data(as(X, "dgCMatrix"), prop, alpha=alpha)
        expect_equal(dense.out, sparse.out1)
        sparse.out2 <- DropletUtils:::.compute_multinom_prob_data(as(X, "dgTMatrix"), prop, alpha=alpha)
        expect_equal(dense.out, sparse.out2)

        # Checking it's the same compared to the reference.
        ref <- apply(X, 2, FUN=ddirmultinom, prob=prop, alpha=alpha)
        combined <- dense.out + DropletUtils:::.compute_multinom_prob_rest(colSums(X), alpha=alpha)
        expect_equal(ref, combined)
    }
})

test_that("p-value calculations are correct without alpha", {
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

test_that("p-value calculations are correct with alpha", {
    # Simulating some Dirichlet-multinomial probabilities. 
    SIMSTUFF <- function(nbarcodes, ngenes, alpha) { 
        totals <- sample(1000, nbarcodes, replace=TRUE)
        ambient <- runif(ngenes)
        true.prob <- ambient * alpha

        probs <- numeric(length(totals))
        for (i in seq_along(totals)) { 
            prob <- rgamma(ngenes, true.prob) # pre-simulation of the Dirichlet
            sim <- rmultinom(1, totals[i], prob=prob)
            probs[i] <- sum(lgamma(sim + true.prob) - lfactorial(sim) - lgamma(true.prob))
        }
        return(list(totals=totals, probs=probs, ambient=ambient))
    }

    # Setting up the reference function.
    REFFUN <- function(totals, probs, ambient, iter, alpha) {
        n.above <- integer(length(totals))
        alpha.ambient <- ambient*alpha

        for (it in seq_len(iter)) {
            test.ambient <- rgamma(length(ambient), alpha.ambient)
            test.ambient <- test.ambient/sum(test.ambient)

            cumamb <- cumsum(test.ambient)
            logamb <- log(test.ambient)
            sumamb <- sum(test.ambient)

            vec <- runif(max(totals), 0, sumamb)
            chosen <- findInterval(vec, cumamb)+1L

            for (i in seq_along(totals)) {
                curchosen <- chosen[seq_len(totals[i])]
                curprof <- tabulate(curchosen, nbins=length(test.ambient))
                curP <- sum(lgamma(alpha.ambient + curprof) - lfactorial(curprof) - lgamma(alpha.ambient))

                if (curP <= probs[i]) { 
                    n.above[i] <- n.above[i] + 1L
                }
            }
        }

        return(n.above)
    }

    # Checking the equivalence of our two implementations.
    set.seed(0)
    sim <- SIMSTUFF(100, 20, alpha=10)
    totals <- sim$totals
    probs <- sim$probs
    ambient.prof <- sim$ambient

    set.seed(100)
    ref <- REFFUN(totals, probs, ambient.prof, iter=100, alpha=10)
    set.seed(100)
    stats <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=100, alpha=10)
    expect_identical(stats, ref)

    # More barcodes (especially checking what happens with two barcodes of the same total).
    set.seed(0)
    sim <- SIMSTUFF(1001, 20, alpha=5)
    totals <- sim$totals
    probs <- sim$probs
    ambient.prof <- sim$ambient

    set.seed(100)
    ref <- REFFUN(totals, probs, ambient.prof, iter=100, alpha=5)
    set.seed(100)
    stats <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=100, alpha=5)
    expect_identical(stats, ref)

    # More genes.
    set.seed(0)
    sim <- SIMSTUFF(100, 200, alpha=2)
    totals <- sim$totals
    probs <- sim$probs
    ambient.prof <- sim$ambient

    set.seed(100)
    ref <- REFFUN(totals, probs, ambient.prof, iter=100, alpha=2)
    set.seed(100)
    stats <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=100, alpha=2)
    expect_identical(stats, ref)

    # More iterations.
    set.seed(0)
    sim <- SIMSTUFF(100, 20, alpha=20)
    totals <- sim$totals
    probs <- sim$probs
    ambient.prof <- sim$ambient

    set.seed(100)
    ref <- REFFUN(totals, probs, ambient.prof, iter=500, alpha=20)
    set.seed(100)
    stats <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=500, alpha=20)
    expect_identical(stats, ref)

    # Checking for uniformity:
#    set.seed(0)
#    sim <- SIMSTUFF(10000, 1000, alpha=5)
#    totals <- sim$totals
#    probs <- sim$probs
#    ambient.prof <- sim$ambient
#    set.seed(100)
#    stats <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=1000, alpha=5)
})

test_that("emptyDrops runs to completion", {
    # Mocking up some counts.
    set.seed(1000)
    my.counts <- DropletUtils:::simCounts()
    limit <- 100
    e.out <- emptyDrops(my.counts, lower=limit, alpha=Inf)
   
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
    e.out2 <- emptyDrops(my.counts, lower=limit, test.ambient=TRUE, alpha=Inf)
    expect_identical(e.out$Total, e.out2$Total)
    expect_identical(e.out$LogProb[totals>limit], e.out2$LogProb[totals>limit])
    expect_true(all(!is.na(e.out2$LogProb[totals>0])))

    # Checking the ignore argument works correctly.
    set.seed(1001)
    e.out3a <- emptyDrops(my.counts, lower=limit, ignore=200, alpha=Inf)
    e.out3a$FDR <- NULL
    set.seed(1001)
    survivors <- totals <= 100 | totals > 200
    new.counts <- my.counts[,survivors]
    e.out3b <- emptyDrops(new.counts, lower=limit, ignore=NULL, alpha=Inf)
    e.out3b$FDR <- NULL
    expect_equal(e.out3a[survivors,], e.out3b) 

    # Checking retention options.
    K <- barcodeRanks(my.counts, lower=limit)$knee
    expect_true(all(e.out$FDR[totals >= K]==0))
    expect_true(!all(e.out$FDR[totals < K]==0))

    e.outK <- emptyDrops(my.counts, retain=K*0.6, alpha=Inf)
    expect_true(all(e.outK$FDR[totals >= K*0.6]==0))
    expect_true(!all(e.outK$FDR[totals < K*0.6]==0))
})

test_that("emptyDrops works correctly with alpha estimation", {
    # Mocking up some counts.
    set.seed(7000)
    my.counts <- DropletUtils:::simCounts() * 2
    limit <- 100
    e.out <- emptyDrops(my.counts, lower=limit, alpha=NULL)

    # Pulling out the alpha.
    totals <- Matrix::colSums(my.counts)
    discarded <- totals <= limit
    lostmat <- my.counts[,discarded,drop=FALSE]
    ambient.prop <- metadata(e.out)$ambient
    alpha0 <- metadata(e.out)$alpha

    alphaT <- DropletUtils:::.estimate_alpha(as(lostmat, "dgTMatrix"), ambient.prop, totals[discarded])
    expect_identical(alpha0, alphaT)

    # Checking that the total probability calculations are correct.
    valid <- which(totals > limit)
    collected <- numeric(nrow(e.out))
    for (v in valid) { 
        collected[v] <- ddirmultinom(my.counts[,v], prob=ambient.prop, alpha=alpha0)
    }
    expect_equal(e.out$LogProb[valid], collected[valid])

    # Checking that it responds to the specified value of alpha.
    e.out2 <- emptyDrops(my.counts, lower=limit, alpha=5)
    expect_false(isTRUE(all.equal(e.out2$PValue, e.out$PValue)))
})
