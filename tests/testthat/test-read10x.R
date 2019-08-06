# Testing the read10xCounts function.
# library(DropletUtils); library(testthat); source("test-read10x.R")

set.seed(2000)
library(Matrix)

# Mocking up some 10X genomics output.
my.counts <- matrix(rpois(1000, lambda=5), ncol=10, nrow=100)
my.counts <- as(my.counts, "dgCMatrix")

ngenes <- nrow(my.counts)
gene.ids <- paste0("GENE", seq_len(ngenes))
gene.symb <- paste0(sample(LETTERS, replace=TRUE, ngenes),
                    sample(LETTERS, replace=TRUE, ngenes),
                    sample(LETTERS, replace=TRUE, ngenes), "-",
                    sample(9, replace=TRUE, ngenes))

cell.ids <- paste0("BARCODE-", seq_len(ncol(my.counts)))

test_that("read10xCounts works correctly for sparse counts, version < 3", {
    tmpdir <- tempfile()
    write10xCounts(path=tmpdir, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids)

    # Reading it in.
    sce10x <- read10xCounts(tmpdir)
    alt.counts <- my.counts
    rownames(alt.counts) <- gene.ids
    colnames(alt.counts) <- NULL

    expect_equal(counts(sce10x), alt.counts)
    expect_identical(rowData(sce10x)$ID, gene.ids)
    expect_identical(rowData(sce10x)$Symbol, gene.symb)
    expect_identical(sce10x$Sample, rep(tmpdir, ncol(my.counts)))
    expect_identical(sce10x$Barcode, cell.ids)

    # Reading it in, twice; and checking it makes sense.
    sce10x2 <- read10xCounts(c(tmpdir, tmpdir))
    ref <- sce10x
    colnames(ref) <- NULL
    ref2 <- cbind(ref, ref)
    int_metadata(ref2) <- int_metadata(ref)
    expect_equal(ref2, sce10x2)

    # Checking that column names work.
    sce10x3 <- read10xCounts(tmpdir, col.names=TRUE)
    expect_identical(colnames(sce10x3), sce10x3$Barcode)
    sce10x4 <- read10xCounts(c(tmpdir, tmpdir), col.names=TRUE)
    expect_identical(colnames(sce10x4), paste0(rep(1:2, each=ncol(sce10x3)), "_", colnames(sce10x3)))
})

test_that("read10xCounts works for sparse counts with odd inputs", {
    # Checking that we are robust to odd symbols in the gene names.
    tmpdir <- tempfile()
    gene.symb2 <- paste0(gene.symb, sample(c("#", "'", '"', ""), length(gene.ids), replace=TRUE))
    write10xCounts(path=tmpdir, my.counts, gene.id=gene.ids, gene.symbol=gene.symb2, barcodes=cell.ids)
    sce10x <- read10xCounts(tmpdir)

    expect_identical(assay(sce10x, withDimnames=FALSE), my.counts)
    expect_identical(colData(sce10x)$Barcode, cell.ids)
    expect_identical(rowData(sce10x)$ID, gene.ids)
    expect_identical(rowData(sce10x)$Symbol, gene.symb2)

    # Checking that we are robust to names in the inputs. 
    tmpdir2 <- tempfile()
    write10xCounts(path=tmpdir2, my.counts, gene.id=gene.ids, gene.symbol=gene.symb2, barcodes=cell.ids)

    sce10x2 <- read10xCounts(c(A=tmpdir, B=tmpdir2))
    expect_identical(assay(sce10x2), cbind(assay(sce10x), assay(sce10x)))

    expect_identical(sce10x2$Barcode, rep(colData(sce10x)$Barcode, 2))
    expect_identical(unname(sce10x2$Sample), rep(c(tmpdir, tmpdir2), each=ncol(sce10x)))
    expect_identical(names(sce10x2$Sample), rep(c("A", "B"), each=ncol(sce10x)))

    expect_identical(rowData(sce10x2)$ID, rowData(sce10x)$ID)
    expect_identical(rowData(sce10x2)$Symbol, gene.symb2)
})

test_that("read10xCounts works correctly for sparse counts, version >= 3", {
    tmpdir <- tempfile()
    gene.type <- sample(LETTERS, ngenes, replace=TRUE)
    write10xCounts(path=tmpdir, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, 
        gene.type=gene.type, barcodes=cell.ids, version="3")

    sce10x <- read10xCounts(tmpdir)
    alt.counts <- my.counts
    rownames(alt.counts) <- gene.ids
    colnames(alt.counts) <- NULL

    expect_equal(counts(sce10x), alt.counts)
    expect_identical(rowData(sce10x)$ID, gene.ids)
    expect_identical(rowData(sce10x)$Symbol, gene.symb)
    expect_identical(rowData(sce10x)$Type, gene.type)
    expect_identical(sce10x$Sample, rep(tmpdir, ncol(my.counts)))
    expect_identical(sce10x$Barcode, cell.ids)
})

test_that("read10xCounts works correctly for HDF5 counts, version < 3", {
    tmph5 <- tempfile(fileext=".h5")
    write10xCounts(path=tmph5, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids)
        
    # Reading it in.
    sce10x <- read10xCounts(tmph5)
    alt.counts <- as.matrix(my.counts)
    dimnames(alt.counts) <- NULL

    expect_s4_class(counts(sce10x, withDimnames=FALSE), "DelayedMatrix")
    expect_equal(as.matrix(counts(sce10x, withDimnames=FALSE)), alt.counts)
    expect_identical(rowData(sce10x)$ID, gene.ids)
    expect_identical(rowData(sce10x)$Symbol, gene.symb)
    expect_identical(sce10x$Sample, rep(tmph5, ncol(my.counts)))
    expect_identical(sce10x$Barcode, cell.ids)

    # Reading it in, twice; and checking it makes sense.
    sce10x2 <- read10xCounts(c(tmph5, tmph5))
    ref <- sce10x
    colnames(ref) <- NULL
    ref <- cbind(ref, ref)
    expect_identical(colData(ref), colData(sce10x2))
    expect_identical(rowData(ref), rowData(sce10x2))
    expect_identical(as.matrix(assay(ref)), as.matrix(assay(sce10x2)))
})

test_that("read10xCounts works correctly for HDF5 counts, version >= 3", {
    tmph5 <- tempfile(fileext=".h5")
    gene.type <- sample(LETTERS, ngenes, replace=TRUE)
    write10xCounts(path=tmph5, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, gene.type=gene.type, barcodes=cell.ids, version="3")
        
    sce10x <- read10xCounts(tmph5)
    alt.counts <- as.matrix(my.counts)
    dimnames(alt.counts) <- NULL

    expect_s4_class(counts(sce10x, withDimnames=FALSE), "DelayedMatrix")
    expect_equal(as.matrix(counts(sce10x, withDimnames=FALSE)), alt.counts)

    expect_identical(rowData(sce10x)$ID, gene.ids)
    expect_identical(rowData(sce10x)$Symbol, gene.symb)
    expect_identical(rowData(sce10x)$Type, gene.type)

    expect_identical(sce10x$Sample, rep(tmph5, ncol(my.counts)))
    expect_identical(sce10x$Barcode, cell.ids)
})

