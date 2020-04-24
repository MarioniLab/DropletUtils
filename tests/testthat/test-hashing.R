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

REF <- function(y, p, pseudo) {
    ncells <- ncol(y)
    nhto <- nrow(y)

    ref.best <- ref.second <- integer(ncells)
    ref.fc <- ref.fc2 <- numeric(ncells)

    for (i in seq_len(ncells)) {
        Y <- y[,i]

        adj <- Y/p
        if (nhto>=5) {
            scaling <- median(adj)
        } else if (nhto==4) {
            scaling <- sort(adj, decreasing=TRUE)[3]
        } else {
            scaling <- min(adj)
        }

        ambient <- scaling * p
        Y0 <- pmax(0, Y - ambient)
        PS <- max(pseudo, mean(ambient))
        Y0 <- Y0 + PS

        o <- order(Y0, decreasing=TRUE)
        ref.best[i] <- o[1]
        ref.second[i] <- o[2]
        
        sorted <- Y0[o]
        ref.fc[i] <- sorted[1]/sorted[2]
        ref.fc2[i] <- sorted[2]/PS
    }

    list(Best=ref.best, Second=ref.second, FC=ref.fc, FC2=ref.fc2)   
}

test_that("hashed_deltas works as expected", {
    # Cycling across pseudo.
    for (PSEUDO in c(1, 3, 5)) {
        p <- runif(nhto)
        output <- DropletUtils:::hashed_deltas(y, p, pseudo=PSEUDO)
        ref <- REF(y, p, PSEUDO)

        expect_identical(output$Best, ref$Best-1L)
        expect_identical(output$Second, ref$Second-1L)
        expect_equal(output$FC, ref$FC)
        expect_equal(output$FC2, ref$FC2)
    }

    # Cycling across different number of samples.
    for (N in 3:6) {
        z <- y[1:N,,drop=FALSE]
        q <- runif(N)
        pseudo <- 1

        output <- DropletUtils:::hashed_deltas(z, q, pseudo=pseudo)
        ref <- REF(z, q, pseudo)

        expect_identical(output$Best, ref$Best-1L)
        expect_identical(output$Second, ref$Second-1L)
        expect_equal(output$FC, ref$FC)
        expect_equal(output$FC2, ref$FC2)
    }
})

test_that("hashed_deltas falls back when there are very few samples", {
    p <- runif(2)
    pseudo <- 2
    output <- DropletUtils:::hashed_deltas(y[1:2,], p, pseudo=pseudo)

    ref <- REF(y[1:2,], p, pseudo)
    expect_equal(output$Best + 1, ref$Best)
    expect_equal(output$FC, ref$FC)

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
    expect_false(any(out$Doublet & out$Confident))
    
    expect_true(min(out$LogFC2[out$Doublet]) > max(out$LogFC2[!out$Doublet]))
    expect_true(min(out$LogFC[out$Confident]) > max(out$LogFC[!out$Doublet & !out$Confident]))

    # Testing against the known truth.
    expect_identical(out$Best, true.sample)
    expect_equal(out$Second[1:ndoub], next.sample[1:ndoub])
    expect_true(min(out$LogFC2[1:ndoub]) > max(out$LogFC2[-(1:ndoub)]))

    # Works with mixture models.
    out <- hashedDrops(y, doublet.mixture=TRUE)
    expect_true(min(out$LogFC2[out$Doublet]) > max(out$LogFC2[!out$Doublet]))
    expect_true(min(out$LogFC[out$Confident]) > max(out$LogFC[!out$Doublet & !out$Confident]))

    # Works with the ambient estimation turned off.
    out <- hashedDrops(y, assume.constant=TRUE)
    expect_true(min(out$LogFC2[out$Doublet]) > max(out$LogFC2[!out$Doublet]))
    expect_true(min(out$LogFC[out$Confident]) > max(out$LogFC[!out$Doublet & !out$Confident]))
})
