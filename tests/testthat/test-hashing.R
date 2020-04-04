# Tests for hashedDrops.
# library(testthat); library(DropletUtils); source("test-hashing.R")

# Mocking up an example dataset with 10 HTOs and 10% doublets.
set.seed(90000)
ncells <- 1000
nhto <- 10
y <- matrix(rpois(ncells*nhto, 50), nrow=nhto)
true.sample <- sample(nhto, ncells, replace=TRUE)
y[cbind(true.sample, seq_len(ncells))] <- 1000

ndoub <- ncells/10
next.sample <- (true.sample[1:ndoub]  + 1) %% nrow(y)
next.sample[next.sample==0] <- nrow(y)
y[cbind(next.sample, seq_len(ndoub))] <- 500

test_that("hashed_deltas works as expected", {
    p <- runif(nhto)
    p <- p/mean(p)
    pseudo <- 1

    output <- DropletUtils:::hashed_deltas(y, p, pseudo=pseudo)

    ref.best <- ref.second <- integer(ncells)
    ref.fc <- ref.fc2 <- numeric(ncells)
    for (i in seq_len(ncells)) {
        current <- y[,i]/p + pseudo

        o <- order(current, decreasing=TRUE)
        ref.best[i] <- o[1]
        ref.second[i] <- o[2]
        
        sorted <- current[o]
        lfc <- diff(log2(sorted))
        ref.fc[i] <- -lfc[1]
        ref.fc2[i] <- -lfc[2]
    }

    expect_identical(output$Best, ref.best-1L)
    expect_identical(output$Second, ref.second-1L)
    expect_equal(output$FC, 2^ref.fc)
    expect_equal(output$FC2, 2^ref.fc2)
})

test_that("hashed_deltas falls back when there are very few samples", {
    p <- runif(2)
    p <- p/mean(p)
    pseudo <- 2
    output <- DropletUtils:::hashed_deltas(y[1:2,], p, pseudo=pseudo)

    adj <- y[1:2,]/p
    expect_equal(output$Best + 1, max.col(t(adj)))
    expect_equal(output$FC, 2^abs( log2(adj[1,] + pseudo) - log2(adj[2,] + pseudo) ))

    expect_true(all(is.na(output$Second)))
    expect_true(all(is.na(output$FC2)))

    # Works with just 1 sample.
    output <- DropletUtils:::hashed_deltas(y[1,,drop=FALSE], p[1], pseudo=pseudo)

    expect_true(all(output$Best==0))
    expect_true(all(is.na(output$FC)))
    expect_true(all(is.na(output$Second)))
    expect_true(all(is.na(output$FC2)))

    # Works with, believe it or not, no samples!
    output <- DropletUtils:::hashed_deltas(y[0,,drop=FALSE], p[0], pseudo=pseudo)

    expect_true(all(is.na(output$Best)))
    expect_true(all(is.na(output$FC)))
    expect_true(all(is.na(output$Second)))
    expect_true(all(is.na(output$FC2)))
})

test_that("hashedDrops works as expected", {
    out <- hashedDrops(y)
    expect_identical(out$Total, colSums(y))
    expect_identical(order(out$NMAD.1to2), order(out$LogFC.1to2))
    expect_identical(order(out$NMAD.2to3), order(out$LogFC.2to3))

    # Testing against the known truth.
    expect_identical(out$Best, true.sample)
    expect_equal(out$Second[1:ndoub], next.sample[1:ndoub])
    expect_true(min(out$LogFC.2to3[1:ndoub]) > max(out$LogFC.2to3[-(1:ndoub)]))
})
