# Testing the read10xCounts function.
# library(DropletUtils); library(testthat); source("test-read10x.R")

set.seed(1000)
library(Matrix)
tmpdir <- tempfile()

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
    
test_that("write10xCounts works correctly", {
    write10xCounts(path=tmpdir, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids)
    expect_identical(sort(list.files(tmpdir)), c("barcodes.tsv", "genes.tsv", "matrix.mtx"))
    all.sizes <- file.info(list.files(tmpdir, full=TRUE))$size

    # Checking overwrite.
    expect_error(write10xCounts(path=tmpdir, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids),
                 "specified 'path' already exists", fixed=TRUE)
    write10xCounts(path=tmpdir, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids, overwrite=TRUE)
    expect_identical(all.sizes, file.info(list.files(tmpdir, full=TRUE))$size)

    # Checking lengths.
    expect_error(write10xCounts(path=tmpdir, my.counts, gene.id=gene.ids, gene.symbol=gene.symb), "barcodes")
    expect_error(write10xCounts(path=tmpdir, my.counts, barcodes=cell.ids, gene.symbol=gene.symb), "lengths of 'gene.id' and 'gene.symbol'")
    expect_error(write10xCounts(path=tmpdir, my.counts, barcodes=cell.ids, gene.id=gene.ids, gene.symbol=""), "lengths of 'gene.id' and 'gene.symbol'")

    all.sizes <- file.info(list.files(tmpdir, full=TRUE))$size # files should still be there after all those errors.
    expect_identical(all.sizes, file.info(list.files(tmpdir, full=TRUE))$size)

    # Checking default arguments.
    new.counts <- my.counts
    rownames(new.counts) <- gene.ids
    colnames(new.counts) <- cell.ids
    write10xCounts(path=tmpdir, new.counts, gene.symbol=gene.symb, overwrite=TRUE)
    expect_identical(all.sizes, file.info(list.files(tmpdir, full=TRUE))$size)
})

test_that("read10xCounts works correctly", {
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
    ref <- BiocGenerics::cbind(ref, ref)
    expect_equal(ref, sce10x2)

    # Checking that column names work.
    sce10x <- read10xCounts(tmpdir, col.names=TRUE)
    expect_identical(colnames(sce10x), sce10x$Barcode)
    sce10x <- read10xCounts(c(tmpdir, tmpdir), col.names=TRUE)
    expect_identical(colnames(sce10x), NULL)
})



