# Testing the downsampleMatrix function.
# library(DropletUtils); library(testthat); source("test-downsample.R")

CHECKFUN <- function(input, prop) {
    out <- downsampleMatrix(input, prop)
    expect_identical(colSums(out), round(colSums(input)*prop))
    expect_true(all(out <= input))
    return(invisible(NULL))
}

CHECKSUM <- function(input, prop) {
    expect_equal(sum(downsampleMatrix(input, prop, bycol=FALSE)), round(prop*sum(input)))
    return(invisible(NULL))
}

test_that("downsampling from a count matrix gives expected sums", {
    # Vanilla run.
    set.seed(0)
    ncells <- 100
    u1 <- matrix(rpois(20000, 5), ncol=ncells)
    u2 <- matrix(rpois(20000, 1), ncol=ncells)

    set.seed(100)
    CHECKFUN(u1, 0.111) # Avoid problems with different rounding of 0.5.
    CHECKFUN(u1, 0.333)
    CHECKFUN(u1, 0.777)
    
    CHECKSUM(u1, 0.111) 
    CHECKSUM(u1, 0.333)
    CHECKSUM(u1, 0.777)

    set.seed(101)
    CHECKFUN(u2, 0.111) # Avoid problems with different rounding of 0.5.
    CHECKFUN(u2, 0.333)
    CHECKFUN(u2, 0.777)

    CHECKSUM(u2, 0.111) 
    CHECKSUM(u2, 0.333)
    CHECKSUM(u2, 0.777)

    # Checking double-precision inputs.
    v1 <- u1
    storage.mode(v1) <- "double"
    set.seed(200)
    CHECKFUN(v1, 0.111)
    CHECKFUN(v1, 0.333)
    CHECKFUN(v1, 0.777)

    CHECKSUM(v1, 0.111)
    CHECKSUM(v1, 0.333)
    CHECKSUM(v1, 0.777)

    v2 <- u2
    storage.mode(v2) <- "double"
    set.seed(202)
    CHECKFUN(v2, 0.111)
    CHECKFUN(v2, 0.333)
    CHECKFUN(v2, 0.777)
    
    CHECKSUM(v2, 0.111)
    CHECKSUM(v2, 0.333)
    CHECKSUM(v2, 0.777)

    # Checking vectors of proportions.
    set.seed(300)
    CHECKFUN(u1, runif(ncells))
    CHECKFUN(u1, runif(ncells, 0, 0.5))
    CHECKFUN(u1, runif(ncells, 0.1, 0.2))

    set.seed(303)
    CHECKFUN(u2, runif(ncells))
    CHECKFUN(u2, runif(ncells, 0, 0.5))
    CHECKFUN(u2, runif(ncells, 0.1, 0.2))

    # Checking sparse matrix inputs.
    library(Matrix)
    w1 <- as(v1, "dgCMatrix")
    set.seed(400)
    CHECKFUN(w1, 0.111)
    CHECKFUN(w1, 0.333)
    CHECKFUN(w1, 0.777)

    CHECKSUM(w1, 0.111)
    CHECKSUM(w1, 0.333)
    CHECKSUM(w1, 0.777)

    w2 <- as(v2, "dgCMatrix")
    set.seed(404)
    CHECKFUN(w2, 0.111)
    CHECKFUN(w2, 0.333)
    CHECKFUN(w2, 0.777)

    CHECKSUM(w2, 0.111)
    CHECKSUM(w2, 0.333)
    CHECKSUM(w2, 0.777)
})

test_that("downsampling from a count matrix gives expected margins", {
    # Checking that the sampling scheme is correct (as much as possible).
    set.seed(500)
    known <- matrix(1:5, nrow=5, ncol=10000)
    prop <- 0.51
    truth <- known[,1]*prop
    out <- downsampleMatrix(known, prop)
    expect_true(all(abs(rowMeans(out)/truth - 1) < 0.1)) # Less than 10% error on the estimated proportions.

    out <- downsampleMatrix(known, prop, bycol=FALSE) # Repeating by column.
    expect_true(all(abs(rowMeans(out)/truth - 1) < 0.1)) 

    # Repeating with larger counts.
    known <- matrix(1:5*100, nrow=5, ncol=10000)
    prop <- 0.51
    truth <- known[,1]*prop
    out <- downsampleMatrix(known, prop)
    expect_true(all(abs(rowMeans(out)/truth - 1) < 0.01)) # Less than 1% error on the estimated proportions.

    out <- downsampleMatrix(known, prop, bycol=FALSE)
    expect_true(all(abs(rowMeans(out)/truth - 1) < 0.01)) 

    # Checking the column sums when bycol=FALSE.
    known <- matrix(100, nrow=1000, ncol=10)
    out <- downsampleMatrix(known, prop, bycol=FALSE)
    expect_true(all(abs(colMeans(out)/colMeans(known)/prop - 1) < 0.01))
})

