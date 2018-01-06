# Testing the downsampleMatrix function.
# library(DropletUtils); library(testthat); source("test-downsample.R")

CHECKFUN <- function(input, prop) {
    out <- downsampleMatrix(input, prop)
    expect_identical(colSums(out), round(colSums(input)*prop))
    expect_true(all(out <= input))
    return(invisible(NULL))
}

CHECKSUM <- function(input, prop) {
    out <- downsampleMatrix(input, prop, bycol=FALSE) 
    expect_equal(sum(out), round(prop*sum(input)))
    expect_true(all(out <= input))
    return(invisible(NULL))
}

test_that("downsampling from a count matrix gives expected sums", {
    # Vanilla run.
    set.seed(0)
    ncells <- 100
    u1 <- matrix(rpois(20000, 5), ncol=ncells)
    u2 <- matrix(rpois(20000, 1), ncol=ncells)

    set.seed(100)
    for (down in c(0.111, 0.333, 0.777)) { # Avoid problems with different rounding of 0.5.
        CHECKFUN(u1, down) 
        CHECKSUM(u1, down) 
    }

    set.seed(101)
    for (down in c(0.111, 0.333, 0.777)) { # Avoid problems with different rounding of 0.5.
        CHECKFUN(u2, down) 
        CHECKSUM(u2, down) 
    }

    # Checking double-precision inputs.
    v1 <- u1
    storage.mode(v1) <- "double"
    set.seed(200)
    for (down in c(0.111, 0.333, 0.777)) { 
        CHECKFUN(v1, down) 
        CHECKSUM(v1, down) 
    }

    v2 <- u2
    storage.mode(v2) <- "double"
    set.seed(202)
    for (down in c(0.111, 0.333, 0.777)) { 
        CHECKFUN(v2, down) 
        CHECKSUM(v2, down) 
    }

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
    for (down in c(0.111, 0.333, 0.777)) { 
        CHECKFUN(w1, down) 
        CHECKSUM(w1, down) 
    }

    w2 <- as(v2, "dgCMatrix")
    set.seed(404)
    for (down in c(0.111, 0.333, 0.777)) { 
        CHECKFUN(w2, down) 
        CHECKSUM(w2, down) 
    }

    # Checking that downsampling preserves relative abundances.
    set.seed(500)
    X <- matrix(1:4*100, ncol=500, nrow=4)
    Y <- downsampleMatrix(X, prop=0.11)
    expect_true(all(abs(rowMeans(Y) - 0.11*rowMeans(X)) < 1))
    Y <- downsampleMatrix(X, prop=0.55)
    expect_true(all(abs(rowMeans(Y) - 0.55*rowMeans(X)) < 1))
    Y <- downsampleMatrix(X, prop=0.11, bycol=FALSE)
    expect_true(all(abs(rowMeans(Y) - 0.11*rowMeans(X)) < 1))
    Y <- downsampleMatrix(X, prop=0.55, bycol=FALSE)
    expect_true(all(abs(rowMeans(Y) - 0.55*rowMeans(X)) < 1))

    # Checking silly inputs.
    expect_equal(downsampleMatrix(u1[0,,drop=FALSE], prop=0.5), u1[0,,drop=FALSE])
    expect_equal(downsampleMatrix(v1[0,,drop=FALSE], prop=0.5), v1[0,,drop=FALSE])
    expect_equal(downsampleMatrix(w1[0,,drop=FALSE], prop=0.5), w1[0,,drop=FALSE])

    expect_equal(downsampleMatrix(u1[,0,drop=FALSE], prop=0.5), u1[,0,drop=FALSE])
    expect_equal(downsampleMatrix(v1[,0,drop=FALSE], prop=0.5), v1[,0,drop=FALSE])
    expect_equal(downsampleMatrix(w1[,0,drop=FALSE], prop=0.5), w1[,0,drop=FALSE])

    expect_equal(downsampleMatrix(u1[0,0,drop=FALSE], prop=0.5), u1[0,0,drop=FALSE])
    expect_equal(downsampleMatrix(v1[0,0,drop=FALSE], prop=0.5), v1[0,0,drop=FALSE])
    expect_equal(downsampleMatrix(w1[0,0,drop=FALSE], prop=0.5), w1[0,0,drop=FALSE])
})

set.seed(500)
test_that("downsampling from a count matrix gives expected margins", {
    # Checking that the sampling scheme is correct (as much as possible).
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

#####################################################
#####################################################
#####################################################

library(Matrix)
set.seed(501)
test_that("downsampling from the reads yields correct results", {
    barcode <- 4L
    tmpdir <- tempfile()
    dir.create(tmpdir)
    out.paths <- DropletUtils:::sim10xMolInfo(tmpdir, nsamples=1, ngenes=100, swap.frac=0, barcode.length=barcode) 
   
    # Creating the full matrix, and checking that it's the same when no downsampling is requested.
    collated <- read10xMolInfo(out.paths, barcode)
    all.cells <- sort(unique(collated$data$cell))
    full.tab <- makeCountMatrix(collated$data$gene, collated$data$cell, all.genes=collated$genes)
    colnames(full.tab) <- paste0(colnames(full.tab), "-1")

    out <- downsampleReads(out.paths, barcode, prop=1)
    expect_equal(out, full.tab)

    # Checking that some downsampling has occurred (hard to check the totals, as UMI counts != read counts).
    for (down in 1:4/11) {
        out <- downsampleReads(out.paths, barcode, prop=down)
        expect_true(all(out <= full.tab))
        expect_false(all(out==full.tab))
    }

    # Making it easier to check the totals, by making all UMIs have a read count of 1.
    out.paths <- DropletUtils:::sim10xMolInfo(tmpdir, nsamples=1, ngenes=100, swap.frac=0, barcode.len=barcode, ave.read=0) 
    full.tab <- downsampleReads(out.paths, barcode, prop=1)
    expect_equal(sum(downsampleReads(out.paths, barcode, prop=0.555)), round(0.555*sum(full.tab))) # Again, avoiding rounding differences.
    expect_equal(sum(downsampleReads(out.paths, barcode, prop=0.111)), round(0.111*sum(full.tab)))
    expect_equal(colSums(downsampleReads(out.paths, barcode, prop=0.555, bycol=TRUE)), round(0.555*colSums(full.tab)))
    expect_equal(colSums(downsampleReads(out.paths, barcode, prop=0.111, bycol=TRUE)), round(0.111*colSums(full.tab)))

    # Manually creating files for comparison to downsampleMatrix - this relies on ordered 'gene' and 'cell', 
    # so that the retention probabilities applied to each count are the same across functions.
    ngenes <- 4
    gene.count <- seq_len(ngenes)*100
    ncells <- 200
    nmolecules <- sum(gene.count)*ncells
    out.file <- file.path(tmpdir, "out.h5")

    library(rhdf5)
    h5 <- h5createFile(out.file)
    h5write(rep(seq_len(ncells), each=sum(gene.count)), out.file, "barcode")
    h5write(seq_len(nmolecules), out.file, "umi")
    h5write(rep(rep(seq_len(ngenes)-1L, gene.count), ncells), out.file, "gene")
    h5write(rep(1, nmolecules), out.file, "gem_group")
    h5write(rep(1, nmolecules), out.file, "reads")
    h5write(array(sprintf("ENSG%i", seq_len(ngenes))), out.file, "gene_ids")

    alt <- read10xMolInfo(out.file)
    X <- makeCountMatrix(alt$data$gene, alt$data$cell, all.genes=alt$genes)
    colnames(X) <- paste0(colnames(X), "-1")
    set.seed(100)
    Z <- downsampleMatrix(X, prop=0.11, bycol=FALSE)
    set.seed(100)
    Y <- downsampleReads(out.file, prop=0.11)
    expect_equal(Y, Z)

    set.seed(100)
    Z <- downsampleMatrix(X, prop=0.55)
    set.seed(100)
    Y <- downsampleReads(out.file, prop=0.55, bycol=TRUE)
    expect_equal(Y, Z)

    # Checking behaviour on silly inputs where there are no reads, or no genes.
    ngenes <- 20L
    out.paths <- DropletUtils:::sim10xMolInfo(tmpdir, nsamples=1, nmolecules=0, swap.frac=0, ngenes=ngenes, barcode.length=barcode) 
    out <- downsampleReads(out.paths, barcode, prop=0.5)
    expect_identical(dim(out), c(ngenes, 0L))

    out.paths <- DropletUtils:::sim10xMolInfo(tmpdir, nsamples=1, nmolecules=0, ngenes=0, swap.frac=0, barcode.length=barcode) 
    out <- downsampleReads(out.paths, barcode, prop=0.5)
    expect_identical(dim(out), c(0L, 0L))
})
