# This tests the removeAmbience function.
# library(testthat); library(DropletUtils); source("test-rmamb.R")

ngenes <- 1000
ncells <- 100
ngroups <- 4

# Making up some data.
ambient <- runif(ngenes, 0, 0.1)
sf <- rep(1:25/25+1, each=ngroups)

means <- matrix(0, nrow=ngenes, ncol=ngroups)
features <- list(1:100, 101:200, 201:300, 301:400)
for (i in seq_len(ngroups)) means[features[[i]],i] <- i * 10

groupings <- rep(1:4, length.out=ncells)
cells <- means[,groupings,drop=FALSE]
y <- t(t(cells + ambient * 100) * sf)

test_that("removeAmbience works correctly in a dummy scenario", {
    removed <- removeAmbience(y, ambient=ambient, groups=groupings, size.factors=sf, features=features)

    accumulated <- numeric(ngroups)
    for (i in seq_len(ngroups)) {
        chosen <- groupings==i
        targets <- removed[,chosen]
        expect_true(all(targets[-features[[i]],]==0))

        residuals <- colSums(targets[features[[i]],])
        expect_false(is.unsorted(residuals))

        ratio <- residuals/sf[chosen]
        cv <- sd(ratio)/mean(ratio)
        expect_true(cv < 0.01) # minimal variation from size factors.

        accumulated[i] <- mean(residuals)
    }

    ratio <- (accumulated/1:4)
    cv <- sd(ratio)/mean(ratio)
    expect_true(cv < 0.01)
})

test_that("removeAmbience works correctly with DF inputs", {
    ref <- removeAmbience(y, ambient=ambient, groups=groupings, size.factors=sf, features=features)
    obs <- removeAmbience(y, ambient=ambient, groups=DataFrame(group=groupings), size.factors=sf, features=features)
    expect_identical(ref, obs)
})

library(DelayedArray)

test_that("removeAmbience works correctly with other input matrices", {
    ref <- removeAmbience(y, ambient=ambient, groups=groupings, size.factors=sf, features=features)

    y2 <- as(y, "dgCMatrix")
    obs <- removeAmbience(y2, ambient=ambient, groups=groupings, size.factors=sf, features=features)
    expect_s4_class(obs, "dgCMatrix")
    expect_equivalent(ref, as.matrix(obs))

    y3 <- DelayedArray(y)
    oldb <- getAutoBlockSize()
    for (blocks in c(1000, 50000, 100000)) {
        setAutoBlockSize(blocks)
        obs <- removeAmbience(y3, ambient=ambient, groups=groupings, size.factors=sf, features=features)
        expect_identical(ref, obs)

        sink <- RealizationSink(dim(ref)) 
        alt <- removeAmbience(y, ambient=ambient, groups=groupings, size.factors=sf, features=features, sink=sink)
        expect_s4_class(alt, "DelayedMatrix")
        expect_identical(ref, as.matrix(alt))
    }
    setAutoBlockSize(oldb)
})
