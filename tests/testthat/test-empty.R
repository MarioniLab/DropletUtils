# This tests that the C++ code for emptyDrops does what it says it does.
# library(DropletUtils); library(testthat); source("setup.R"); source("test-empty.R")

test_that("estimateAmbience works correctly", {
    # Mocking up some counts.
    set.seed(1000)
    my.counts <- DropletUtils:::simCounts()
    limit <- 100

    out <- estimateAmbience(my.counts)
    keep <- colSums(my.counts) <= limit
    naive <- rowSums(my.counts[,keep])
    expect_false(is.unsorted(out[order(naive)]))
    expect_identical(naive, estimateAmbience(my.counts, good.turing=FALSE))

    my.counts2 <- my.counts
    rownames(my.counts2) <- NULL
    out2 <- estimateAmbience(my.counts2)
    expect_identical(unname(out), out2)

    # Performs correctly with the rank-based method.
    out.r <- estimateAmbience(my.counts, by.rank=1000)
    expect_false(identical(out, out.r))
    out.r2 <- estimateAmbience(my.counts, by.rank=1000, good.turing=FALSE)
    expect_identical(out.r2, rowSums(my.counts[,rank(-colSums(my.counts), ties.method="first") > 1000]))

    # Deals with zeroes properly.
    my.counts3 <- rbind(0,0,0,my.counts)
    out3 <- estimateAmbience(my.counts3)
    expect_identical(out, tail(out3, -3))
})

test_that("Good-Turing protection works correctly", {
    # Protection gives the same estimates for zeroes when a count
    # of 1 is redistributed somewhere else.
    out <- DropletUtils:::.safe_good_turing(c(0,0,2,4))
    ref <- DropletUtils:::.safe_good_turing(c(0,0,1,2,3))

    expect_equal(sum(ref), 1)
    expect_equal(sum(out), 1)
    expect_false(any(ref==0))
    expect_false(any(out==0))

    expect_equal(out[1:2], ref[1:2])

    # Continues to work with loads of observations.
    out <- DropletUtils:::.safe_good_turing(c(0,0,2,4,10000))
    ref <- DropletUtils:::.safe_good_turing(c(0,0,1,2,3,10000))
    expect_equal(out[1:2], ref[1:2])
})

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
    sparse.out2 <- DropletUtils:::.compute_multinom_prob_data(as(X, "TsparseMatrix"), prop)
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
        sparse.out2 <- DropletUtils:::.compute_multinom_prob_data(as(X, "TsparseMatrix"), prop, alpha=alpha)
        expect_equal(dense.out, sparse.out2)

        # Checking it's the same compared to the reference.
        ref <- apply(X, 2, FUN=ddirmultinom, prob=prop, alpha=alpha)
        combined <- dense.out + DropletUtils:::.compute_multinom_prob_rest(colSums(X), alpha=alpha)
        expect_equal(ref, combined)
    }
})

rdirmultinom <- function(ambient, totals, alpha) {
    true.prob <- ambient * alpha
    collected <- vector("list", length(totals))
    for (i in seq_along(collected)) {
        prob <- rgamma(length(true.prob), true.prob) # pre-simulation of the Dirichlet
        sim <- rmultinom(1, totals[i], prob=prob)
        collected[[i]] <- sim
    }
    collected
}

set.seed(100010010)
test_that("dispersion calculation is close to correct", {
    totals <- sample(10000, 100, replace=TRUE)
    ambient <- runif(2000)
    ambient <- ambient/sum(ambient)
    
    simulated <- rdirmultinom(ambient, totals, alpha=5)
    mat <- do.call(cbind, simulated)
    x <- DropletUtils:::.estimate_alpha(mat, ambient, totals)
    expect_true(abs(x - 5) <= 0.2)

    simulated <- rdirmultinom(ambient, totals, alpha=20)
    mat <- do.call(cbind, simulated)
    x <- DropletUtils:::.estimate_alpha(mat, ambient, totals)
    expect_true(abs(x - 20) <= 1)

    simulated <- rdirmultinom(ambient, totals, alpha=100)
    mat <- do.call(cbind, simulated)
    x <- DropletUtils:::.estimate_alpha(mat, ambient, totals)
    expect_true(abs(x - 100) <= 5)
})

UNICHECKER <- function(values, tol=1.25) {
    for (thresh in c(0.01, 0.05, 0.1, 0.2, 0.5)) { 
        # Generous boundaries for thresholds.
        observed <- mean(values < thresh)
        expect_true(observed < thresh * tol)
        expect_true(observed > thresh / tol)
    }
    invisible(NULL)
}

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

    # Checking for uniformity:
    set.seed(0)
    sim <- SIMSTUFF(10000, 1000)
    totals <- sim$totals
    probs <- sim$probs
    ambient.prof <- sim$ambient

    set.seed(100)
    N <- 10000
    stats <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=N)
    UNICHECKER(stats/N)

    set.seed(100)
    stats2 <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, BPPARAM=safeBPPARAM, iter=N) 
    expect_identical(stats, stats2)

    # Different result for different seeds.
    set.seed(101)
    stats3 <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=N) 
    expect_false(identical(stats, stats3))
    UNICHECKER(stats3/N)

    # Worker splitting behaves for near-zero or no jobs.
    almost_none <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, BPPARAM=safeBPPARAM, iter=1L) 
    expect_true(all(almost_none %in% 0:1))
    none <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, BPPARAM=safeBPPARAM, iter=0L) 
    expect_identical(none, integer(length(totals)))
})

test_that("p-value calculations are correct with alpha", {
    # Simulating some Dirichlet-multinomial probabilities. 
    SIMSTUFF <- function(nbarcodes, ngenes, alpha) { 
        totals <- sample(1000, nbarcodes, replace=TRUE)
        ambient <- runif(ngenes)
        simulated <- rdirmultinom(ambient, totals, alpha)

        true.prob <- ambient * alpha
        probs <- numeric(length(totals))
        for (i in seq_along(simulated)) { 
            sim <- simulated[[i]]
            probs[i] <- sum(lgamma(sim + true.prob) - lfactorial(sim) - lgamma(true.prob))
        }

        list(totals=totals, probs=probs, ambient=ambient)
    }

    # Checking for uniformity.
    set.seed(0)
    sim <- SIMSTUFF(10000, 800, alpha=10)
    totals <- sim$totals
    probs <- sim$probs
    ambient.prof <- sim$ambient

    set.seed(100)
    N <- 10000
    stats <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=N, alpha=10)
    UNICHECKER(stats/N)

    # Same result for multiple cores (set BPPARAM first as it redefines the seed!)
    set.seed(100)
    stats2 <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, BPPARAM=safeBPPARAM, iter=N, alpha=10) 
    expect_identical(stats, stats2)

    # Different result for different seeds.
    set.seed(101)
    stats3 <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=N, alpha=10) 
    expect_false(identical(stats, stats3))
    UNICHECKER(stats3/N)
 
    # Trying with a different alpha.
    set.seed(0)
    sim <- SIMSTUFF(10000, 800, alpha=5)
    totals <- sim$totals
    probs <- sim$probs
    ambient.prof <- sim$ambient

    set.seed(100)
    stats4 <- DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=N, alpha=5)
    UNICHECKER(stats4/N)
    expect_false(identical(stats, stats4))
    
    # Triggering a check on 'alpha'.
    expect_error(DropletUtils:::.permute_counter(totals, probs, ambient.prof, iter=N, alpha=0), "must be positive") 
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
    K <- metadata(barcodeRanks(my.counts, lower=limit))$knee
    expect_true(all(e.out$FDR[totals >= K]==0))
    expect_true(!all(e.out$FDR[totals < K]==0))

    e.outK <- emptyDrops(my.counts, retain=K*0.6, alpha=Inf)
    expect_true(all(e.outK$FDR[totals >= K*0.6]==0))
    expect_true(!all(e.outK$FDR[totals < K*0.6]==0))
})

test_that("emptyDrops ambient tests work correctly", {
    set.seed(1000)
    my.counts <- DropletUtils:::simCounts()
    limit <- 100

    set.seed(1001)
    e.out <- emptyDrops(my.counts, lower=limit, alpha=Inf)
   
    set.seed(1001)
    e.out2 <- emptyDrops(my.counts, lower=limit, test.ambient=TRUE, alpha=Inf)
    expect_identical(e.out$Total, e.out2$Total)

    keep <- e.out$Total > limit
    expect_identical(e.out$LogProb[keep], e.out2$LogProb[keep])
    expect_identical(e.out$PValue[keep], e.out2$PValue[keep])
    expect_identical(e.out$FDR, e.out2$FDR)
    expect_true(all(!is.na(e.out2$LogProb[e.out$Total>0])))

    # Respects the NA setting.
    set.seed(1001)
    e.out3 <- emptyDrops(my.counts, lower=limit, test.ambient=NA, alpha=Inf)
    expect_identical(e.out2$LogProb, e.out3$LogProb)
    expect_identical(e.out2$PValue, e.out3$PValue)
    expect_false(identical(e.out2$FDR, e.out3$FDR))
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

    alphaT <- DropletUtils:::.estimate_alpha(as(lostmat, "TsparseMatrix"), ambient.prop, totals[discarded])
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

test_that("emptyDrops automatically rounds non-integer values", {
    # Mocking up some counts.
    set.seed(7001)
    my.counts <- DropletUtils:::simCounts() * 1.2
    expect_identical(round(my.counts), DropletUtils:::.rounded_to_integer(my.counts))

    set.seed(7002)
    out <- emptyDrops(my.counts)
    set.seed(7002)
    ref <- emptyDrops(round(my.counts))
    expect_identical(out,ref)
})

set.seed(80001)
test_that("emptyDrops works with by.rank=TRUE", {
    my.counts <- DropletUtils:::simCounts() 

    set.seed(999)
    out <- emptyDrops(my.counts, lower=NA, by.rank=1000) # forcing the old lower to not be used.
    expect_true(is.finite(metadata(out)$lower))
    expect_false(metadata(out)$lower==100)
    expect_true(all(is.na(out$LogProb[out$Totals < metadata(out)$lower])))

    set.seed(999)
    ref <- emptyDrops(my.counts, lower=metadata(out)$lower, by.rank=NULL)
    expect_identical(out, ref)

    # Interacts properly with test.ambient=
    set.seed(1000)
    e.out <- emptyDrops(my.counts, lower=NA, by.rank=1000, alpha=Inf)

    set.seed(1000)
    e.out2 <- emptyDrops(my.counts, lower=NA, by.rank=1000, test.ambient=TRUE, alpha=Inf)
    keep <- e.out$Total > metadata(e.out)$lower
    expect_identical(e.out$LogProb[keep], e.out2$LogProb[keep])
    expect_identical(e.out$PValue[keep], e.out2$PValue[keep])
    expect_identical(e.out$FDR, e.out2$FDR)
})

set.seed(80002)
test_that("emptyDrops realizes HDF5Arrays properly", {
    my.counts <- DropletUtils:::simCounts() 

    set.seed(999)
    out <- emptyDrops(my.counts)

    set.seed(999)
    library(HDF5Array)
    da <- as(my.counts, "HDF5Array")
    ref <- emptyDrops(da)
    expect_identical(out, ref)
})

test_that("emptyDrops fails when you don't give it low counts", {
    y <- matrix(rpois(100000, lambda=100), ncol=10000)
    expect_error(emptyDrops(y), "no counts available")
})

set.seed(80003)
test_that("emptyDropsCellRanger does something sensible", {
    my.counts <- DropletUtils:::simCounts(nempty=100000, nlarge=2000, nsmall=1000)
    out <- emptyDropsCellRanger(my.counts)

    s <- colSums(my.counts)
    expect_true(all(is.na(out$FDR[rank(-s) > 3000]))) # known ambients always ignored
    expect_true(all(out$FDR[rank(-s) < 2000] == 0)) # large cells always detected
})
