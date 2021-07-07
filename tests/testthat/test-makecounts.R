# This tests the behaviour of the makeCountMatrix function.
# library(DropletUtils); library(testthat); source("test-makecounts.R")

all.genes <- LETTERS
ngenes <- length(LETTERS)
all.cells <- 1:10
ncells <- length(all.cells) 

test_that("makeCountMatrix works with all input combinations", {
    for (N in c(0, 10, 100, 1000)) { 

    for (G in 1:4) {
        if (G<=2L) {
            genes <- sample(ngenes, N, replace=TRUE)
        } else {
            genes <- sample(all.genes, N, replace=TRUE)
        } 
        if (G%%2==1L) {
            ref.genes <- all.genes
        } else {
            ref.genes <- NULL
        }

    for (C in 1:4) {
        if (C<=2L) {
            cells <- sample(ncells, N, replace=TRUE)
        } else {
            cells <- sample(all.cells, N, replace=TRUE)
        } 
        if (C%%2==1L) {
            ref.cells <- all.cells
        } else {
            ref.cells <- NULL
        }
    
    for (V in 1:2) {
        if (V==1L) {
            vals <- NULL
        } else {
            vals <- rpois(N, lambda=10)
        }
    
    ######################### TESTING BEGINS HERE.

    obs <- DropletUtils:::makeCountMatrix(genes, cells, all.genes=ref.genes, all.cells=ref.cells, value=vals)
    ref <- as.matrix(obs)
    for (i in seq_len(N)) {
        if (is.null(vals)) {
            loss <- 1
        } else {
            loss <- vals[i]
        } 
        ref[genes[i], cells[i]] <- ref[genes[i], cells[i]] - loss
    }
    expect_true(all(ref==0))

    # Can only check dimensions if we've specified them beforehand.    
    if (!is.null(ref.genes)) { 
        expect_identical(nrow(obs), ngenes)
    } 
    if (!is.null(ref.cells)) { 
        expect_identical(ncol(obs), ncells)
    }

    ######################### END TESTING.    

    } # V loop
    } # C loop
    } # G loop
    } # N loop
})

test_that("makeCountMatrix behaves on silly inputs", {
    expect_error(DropletUtils:::makeCountMatrix(1, 1:2), "identical")
    expect_error(DropletUtils:::makeCountMatrix(1, 1, value=1:2), "identical")

    obs <- DropletUtils:::makeCountMatrix(integer(0), integer(0))
    expect_identical(dim(obs), c(0L, 0L))
    obs <- DropletUtils:::makeCountMatrix(integer(0), integer(0), all.genes=all.genes)
    expect_identical(dim(obs), c(ngenes, 0L))
    obs <- DropletUtils:::makeCountMatrix(integer(0), integer(0), all.cells=all.cells)
    expect_identical(dim(obs), c(0L, ncells))

    expect_error(DropletUtils:::makeCountMatrix("a", 1, all.genes=LETTERS), "entry of")
    expect_error(DropletUtils:::makeCountMatrix(100, 1, all.genes=LETTERS), "length of")
    expect_error(DropletUtils:::makeCountMatrix(0, 1), "positive") 

    expect_error(DropletUtils:::makeCountMatrix(1, "a", all.cells=LETTERS), "entry of")
    expect_error(DropletUtils:::makeCountMatrix(1, 100, all.cells=LETTERS), "length of")
    expect_error(DropletUtils:::makeCountMatrix(1, 0), "positive") 
})
