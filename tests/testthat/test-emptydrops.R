# This tests that the C++ code for emptyDrops does what it says it does.
# library(DropletUtils); library(testthat); source("test-emptydrops.R")

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
    ref <- DropletUtils:::.simulate_dev(sim.totals, ambient.prop) 
    expect_equal(sim.LR, ref)             
})

test_that("p-value calculations are correct", {
    ncells <- 10000
    obs.totals <- sample(1000, ncells, replace=TRUE) + 10
    obs.spread <- runif(ncells)
    sim.totals <- seq(0, 2000, by=0.5)
    sim.spread <- runif(length(sim.totals))

    span <- sqrt(2)
    p <- numeric(length(obs.spread))
    limited <- logical(length(obs.spread))
    for (x in seq_along(p)) {
        current <- obs.totals[x]/span <= sim.totals & obs.totals[x]*span >= sim.totals
        nexceed <- sum(sim.spread[current] >= obs.spread[x])
        limited[x] <- nexceed==0L
        p[x] <- (nexceed + 1)/(sum(current)+1)
    }
    
    stats <- DropletUtils:::.compute_P(obs.totals, obs.spread, sim.totals, sim.spread, span=log2(span))
    expect_equal(stats$p.value, p)
    expect_identical(stats$limited, limited)
})

test_that("emptyDrops runs to completion", {
    # Mocking up some counts.
    source(system.file("scripts", "mock_empty.R", package="DropletUtils"))
    limit <- 100
    e.out <- emptyDrops(my.counts, lower=limit, scale=1)
    
    totals <- colSums(my.counts)
    expect_identical(totals, e.out$Total)
    expect_true(all(is.na(e.out$Deviance[totals<=limit])))
    expect_true(all(!is.na(e.out$Deviance[totals>limit])))
    
    K <- findKneePoint(my.counts, lower=limit)
    expect_true(all(e.out$FDR[totals >= K]==0))

    # Checking ambient tests.
    e.out2 <- emptyDrops(my.counts, test.ambient=TRUE)
    expect_identical(e.out$Total, e.out2$Total)
    expect_identical(e.out$Deviance[totals>limit], e.out2$Deviance[totals>limit])
    expect_true(all(!is.na(e.out2$Deviance[totals>0])))

    # Checking scaling options.
    e.out <- emptyDrops(my.counts, scale=0.6)
    expect_true(all(e.out$FDR[totals >= K*0.6]==0))
})
