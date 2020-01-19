# Testing the read10xCounts function.
# library(DropletUtils); library(testthat); source("test-write10x.R")

set.seed(1000)
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

test_that("write10xCounts works correctly for sparse counts, version < 3", {
    tmpdir <- tempfile()
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

test_that("write10xCounts works correctly for sparse counts, version >= 3", {
    tmpdir <- tempfile()
    write10xCounts(path=tmpdir, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids, version="3")
    expect_identical(sort(list.files(tmpdir)), c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"))
    out <- read.table(file.path(tmpdir, "features.tsv.gz"), stringsAsFactors=FALSE, sep="\t")
    expect_identical(out[,3], rep("Gene Expression", ngenes))

    types <- sample(c("Gene Expression", "Antibody", "CUSTOM"), ngenes, replace=TRUE)
    write10xCounts(path=tmpdir, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, gene.type=types, barcodes=cell.ids, version="3", overwrite=TRUE)
    out <- read.table(file.path(tmpdir, "features.tsv.gz"), stringsAsFactors=FALSE, sep="\t")
    expect_identical(out[,3], types)
})

test_that("write10xCounts works correctly for HDF5 counts, version < 3", {
    tmph5 <- tempfile(fileext=".h5")
    write10xCounts(path=tmph5, genome="mm9", my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids)
    all_fields <- rhdf5::h5ls(tmph5)
    expect_identical(all_fields$name, c("mm9", "barcodes", "data", "gene_names", "genes", "indices", "indptr", "shape"))
    all.sizes <- file.info(list.files(tmph5, full=TRUE))$size

    # Checking overwrite.
    expect_error(write10xCounts(path=tmph5, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids),
                 "specified 'path' already exists", fixed=TRUE)
    write10xCounts(path=tmph5, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids, overwrite=TRUE)
    expect_identical(all.sizes, file.info(list.files(tmph5, full=TRUE))$size)
})

test_that("write10xCounts works correctly for HDF5 counts, version >= 3", {
    tmph5 <- tempfile(fileext=".h5")
    write10xCounts(path=tmph5, genome="mm9", my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids, version="3")

    all_fields <- rhdf5::h5ls(tmph5)
    expect_identical(all_fields$name, 
        c("matrix", "barcodes", "data", "features", 
            "_all_tag_keys", "feature_type", "genome", "id", "name", "indices", "indptr", "shape"))

    expect_identical(as.vector(rhdf5::h5read(tmph5, "matrix/features/feature_type")), rep("Gene Expression", ngenes))
    expect_identical(as.vector(rhdf5::h5read(tmph5, "matrix/features/genome")), rep("mm9", ngenes))

    # Overwriting with different genomes and types.
    genomes <- sample(c("mm9", "hg19"), ngenes, replace=TRUE)
    types <- sample(c("Gene Expression", "Antibody", "CUSTOM"), ngenes, replace=TRUE)
    write10xCounts(path=tmph5, genome=genomes, my.counts, 
        gene.id=gene.ids, gene.symbol=gene.symb, gene.type=types, barcodes=cell.ids, version="3", overwrite=TRUE)

    expect_identical(as.vector(rhdf5::h5read(tmph5, "matrix/features/feature_type")), types)
    expect_identical(as.vector(rhdf5::h5read(tmph5, "matrix/features/genome")), genomes)
})
