# This tests that the C++ code does what it says it does.
# library(EmptyDrops); library(testthat); source("test-drops.R")

test_that("simulated deviance calculations are correct", {
    set.seed(100)
    ambient.prop <- runif(1000)
    ambient.prop <- ambient.prop/sum(ambient.prop)

    set.seed(200)
    sim.totals <- seq(1, 2000, by=0.5)
    sim.LR <- numeric(length(sim.totals))
    for (x in seq_along(sim.LR)) {
        cur.means <- ambient.prop*sim.totals[x]
        current <- rpois(length(ambient.prop), lambda=cur.means)
        sim.LR[x] <- sum(current * log(current/cur.means), na.rm=TRUE) + sum(cur.means - current)
    }

    set.seed(200)
    ref <- EmptyDrops:::.simulate_dev(sim.totals, ambient.prop) 
    expect_equal(sim.LR, ref)             
})

test_that("p-value calculations are correct", {
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
