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
    ref.ps <- numeric(ncells)

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
        ref.ps[i] <- PS
    }

    list(Best=ref.best, Second=ref.second, FC=ref.fc, FC2=ref.fc2, Pseudo=ref.ps)  
}

test_that("hashed_deltas works as expected", {
    # Cycling across pseudo.
    for (PSEUDO in c(1, 3, 5)) {
        p <- runif(nhto)
        output <- DropletUtils:::hashed_deltas(y, p, pseudo=PSEUDO, n_expected=1)
        ref <- REF(y, p, PSEUDO)

        expect_identical(drop(output$Best), ref$Best-1L)
        expect_identical(output$Second, ref$Second-1L)
        expect_equal(output$FC, ref$FC)
        expect_equal(output$FC2, ref$FC2)
    }

    # Cycling across different number of samples.
    for (N in 3:6) {
        z <- y[1:N,,drop=FALSE]
        q <- runif(N)
        pseudo <- 1

        output <- DropletUtils:::hashed_deltas(z, q, pseudo=pseudo, n_expected=1)
        ref <- REF(z, q, pseudo)

        expect_identical(drop(output$Best), ref$Best-1L)
        expect_identical(output$Second, ref$Second-1L)
        expect_equal(output$FC, ref$FC)
        expect_equal(output$FC2, ref$FC2)
    }
})

test_that("hashed_deltas falls back when there are very few samples", {
    p <- runif(2)
    pseudo <- 2
    output <- DropletUtils:::hashed_deltas(y[1:2,], p, pseudo=pseudo, n_expected=1)

    ref <- REF(y[1:2,], p, pseudo)
    expect_equal(drop(output$Best) + 1, ref$Best)
    expect_equal(output$FC, ref$FC)

    expect_true(all(is.na(output$Second)))
    expect_true(all(is.na(output$FC2)))

    # Works with just 1 sample.
    output <- DropletUtils:::hashed_deltas(y[1,,drop=FALSE], p[1], pseudo=pseudo, n_expected=1)

    expect_true(all(output$Best==0))
    expect_true(all(is.na(output$FC)))
    expect_true(all(is.na(output$Second)))
    expect_true(all(is.na(output$FC2)))

    # Works with, believe it or not, no samples!
    output <- DropletUtils:::hashed_deltas(y[0,,drop=FALSE], p[0], pseudo=pseudo, n_expected=1)

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

    expect_identical(out$Best, true.sample)
    expect_equal(out$Second[1:ndoub], next.sample[1:ndoub])
    expect_true(min(out$LogFC2[1:ndoub]) > max(out$LogFC2[-(1:ndoub)]))

    # Works with the ambient estimation turned off.
    out <- hashedDrops(y, ambient=rep(1, nrow(y)))
    expect_true(all(metadata(out)$ambient==1))

    expect_true(min(out$LogFC2[out$Doublet]) > max(out$LogFC2[!out$Doublet]))
    expect_true(min(out$LogFC[out$Confident]) > max(out$LogFC[!out$Doublet & !out$Confident]))

    expect_identical(out$Best, true.sample)
    expect_equal(out$Second[1:ndoub], next.sample[1:ndoub])
    expect_true(min(out$LogFC2[1:ndoub]) > max(out$LogFC2[-(1:ndoub)]))
})

test_that("hashedDrops handles low number of tags gracefully", {
    expect_warning(out <- hashedDrops(y[1:2,]), "consider")
    expect_true(all(!is.na(out$LogFC)))
    expect_true(all(!is.na(out$Confident)))
    expect_true(all(is.na(out$LogFC2)))
    expect_true(all(is.na(out$Doublet)))

    out <- hashedDrops(y[1,,drop=FALSE])
    expect_true(all(is.na(out$LogFC)))
    expect_true(all(is.na(out$Confident)))
    expect_true(all(is.na(out$LogFC2)))
    expect_true(all(is.na(out$Doublet)))
})

test_that("hashedDrops works correctly with combinatorial barcodes", {
    mat <- matrix(5, 10, 9)
    mat[c(1, 2, 3), 1] <- sample(c(50, 60, 70))
    mat[c(2, 4, 7), 2] <- sample(c(80, 50, 75))
    mat[c(3, 8, 9), 3] <- sample(c(80, 50, 75))
    mat[c(5, 6), 4] <- c(100, 80)
    mat[c(3, 8), 5] <- c(90, 50)
    mat[7, 6] <- 100
    mat[9, 7] <- 90
    mat[c(1, 3, 5, 7), 8] <- 100

    out <- hashedDrops(mat, combinations=rbind(1:3, c(2,4,7)), ambient=rep(1, nrow(mat)))

    expect_identical(out$Best, c(1:2, rep(NA_integer_, ncol(mat)-2)))
    expect_identical(out$LogFC, rep(c(log2(50/5), 0), c(3, ncol(mat)-3)))

    expect_null(out$Second)

    expect_true(out$Doublet[8])
    expect_true(all(!out$Doublet[-8]))
    expect_identical(out$LogFC2[8], log2(100/5))
    expect_true(all(out$LogFC2[-8]==0))
    
    expect_true(all(out$Confident[1:3]))
    expect_true(all(!out$Confident[-(1:3)]))

    # Same results with unsorted barcodes.
    out2 <- hashedDrops(mat, combinations=rbind(3:1, c(7,2,4)), ambient=rep(1, nrow(mat)))
    expect_identical(out, out2)
})

test_that("edge cases are handled correctly with combinatorial barcodes", {
    # Doublet statistics nullified with insufficient HTOs.
    mat <- matrix(10, 4, 1)
    mat[2:4,1] <- 100

    expect_warning(out <- hashedDrops(mat, combinations=rbind(4:2), ambient=rep(1, nrow(mat))), "consider")
    expect_identical(out$LogFC, log2(100/10))
    expect_identical(out$Best, 1L)
    expect_identical(out$LogFC2, NA_real_)
    expect_identical(out$Doublet, NA)

    # Doublet statistics come back online with just enough HTOs.
    mat <- matrix(10, 7, 1)
    mat[2:4,1] <- 100

    out <- hashedDrops(mat, combinations=rbind(4:2), ambient=rep(1, nrow(mat)))
    expect_identical(out$LogFC, log2(100/10))
    expect_identical(out$Best, 1L)
    expect_identical(out$LogFC2, 0)
    expect_true(!is.na(out$Doublet))

    # What is used to compute the ambient profile?
    library(DropletUtils)
    mat <- cbind(1:20) * 10
    ambient <- 1 + 1:20/1000

    .compute_expected_lfc <- function(SCALING, PSEUDO) {
        log2((180 - SCALING * 1.018  + PSEUDO)/
            (170 - SCALING * 1.017 + PSEUDO))
    }

    for (counter in 17:15) { # Using the last one...
        keep <- 20:counter
        expect_warning(out <- hashedDrops(mat[keep,,drop=FALSE], combinations=rbind(4:2), ambient=ambient[keep]), "consider")
        SCALING <- mat[counter]/ambient[counter]
        PSEUDO <- mean(ambient[keep]) * SCALING
        expect_identical(out$LogFC, .compute_expected_lfc(SCALING, PSEUDO))
    }

    for (counter in 14:9) { # Using the first past the 2*n_expected...
        keep <- 20:counter
        out <- hashedDrops(mat[keep,,drop=FALSE], combinations=rbind(4:2), ambient=ambient[keep])
        SCALING <- mat[14]/ambient[14]
        PSEUDO <- mean(ambient[keep]) * SCALING
        expect_identical(out$LogFC, .compute_expected_lfc(SCALING, PSEUDO))
    }

    for (counter in 8:1) { # an actual median for the rest.
        keep <- 20:counter
        out <- hashedDrops(mat[keep,,drop=FALSE], combinations=rbind(4:2), ambient=ambient[keep])
        SCALING <- median(mat[keep]/ambient[keep])
        PSEUDO <- mean(ambient[keep]) * SCALING
        expect_identical(out$LogFC, .compute_expected_lfc(SCALING, PSEUDO))
    }
})

test_that("hashedDrops works with constant ambience", {
    ref <- hashedDrops(y)
    out <- hashedDrops(y, constant.ambient=TRUE)
    same.col <- c("Total", "Best", "Second", "LogFC")
    expect_identical(ref[,same.col], out[,same.col])

    ref.calc <- REF(y, metadata(ref)$ambient, 5)
    lfc2 <- log2(ref.calc$FC2 * ref.calc$Pseudo / median(ref.calc$Pseudo))
    expect_equal(lfc2, out$LogFC2)

    expect_identical(out$Best, true.sample)
    expect_equal(out$Second[1:ndoub], next.sample[1:ndoub])
    expect_true(min(out$LogFC2[1:ndoub]) > max(out$LogFC2[-(1:ndoub)]))
})
