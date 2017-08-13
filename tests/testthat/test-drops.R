# This tests that the C++ code does what it says it does.
# library(EmptyDrops); library(testthat)

test_that("pvalue calculations are correct", {
    ncells <- 10000
    obs.totals <- sample(1000, ncells, replace=TRUE) + 10
    obs.spread <- runif(ncells)
    sim.totals <- seq(0, 2000, by=0.5)
    sim.spread <- runif(length(sim.totals))

    p <- numeric(length(obs.spread))
    limited <- logical(length(obs.spread))
    fold.tol <- 2^-0.5
    for (x in seq_along(p)) {
        current <- obs.totals[x]*fold.tol <= sim.totals & obs.totals[x]/fold.tol >= sim.totals
        nexceed <- sum(sim.spread[current] >= obs.spread[x])
        limited[x] <- nexceed==0L
        p[x] <- (nexceed + 1)/(sum(current)+1)
    }
    
    stats <- EmptyDrops:::.compute_P(obs.totals, obs.spread, sim.totals, sim.spread, tol=0.5)
    expect_equal(stats$p.value, p)
    expect_identical(stats$limited, limited)
})
