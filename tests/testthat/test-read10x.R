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

    # Adding another dataset with slightly different counts.
    tmpdir2 <- tempfile()
    write10xCounts(path=tmpdir2, my.counts*2, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids)

    sce10x2 <- read10xCounts(tmpdir2)
    expect_identical(assay(sce10x)*2, assay(sce10x2))

    ref <- cbind(sce10x, sce10x2)
    colnames(ref) <- NULL
    combined <- read10xCounts(c(tmpdir, tmpdir2))

    expect_equal(rowData(ref), rowData(combined))
    expect_equal(colData(ref), colData(combined))
    expect_equal(assay(ref), assay(combined))
})

test_that("read10xCounts works correctly with chromosomal positions in the features", {
    tmpdir <- tempfile()
    write10xCounts(path=tmpdir, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids)

    feats <- read.delim(file.path(tmpdir, "genes.tsv"), header=FALSE)
    feats$Type <- sample(c("protein", "foo"), nrow(feats), replace=TRUE)
    feats$Chr <- sample(c("chrA", "chrB", "chrC"), nrow(feats), replace=TRUE)
    feats$Start <- sample(1000, nrow(feats), replace=TRUE)
    feats$End <- feats$Start + sample(100, nrow(feats), replace=TRUE)
    write.table(file=file.path(tmpdir, "genes.tsv"), feats, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

    sce10x <- read10xCounts(tmpdir)
    expect_s4_class(rowRanges(sce10x), "GRanges")
    expect_identical(start(rowRanges(sce10x)), feats$Start)
    expect_identical(rowRanges(sce10x)$Type, feats$Type)

    # Empty seqnames for mitochondria are replaced with chrM.
    replace <- 1:10
    feats$Chr[replace] <- ""
    feats[,2] <- paste0("MT-", feats[,2])
    write.table(file=file.path(tmpdir, "genes.tsv"), feats, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

    sce10x <- read10xCounts(tmpdir)
    expect_identical(unique(as.character(seqnames(rowRanges(sce10x))[1:10])), "chrM")
})

test_that("read10xCounts works correctly for names", {
    tmpdir <- tempfile()
    write10xCounts(path=tmpdir, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids)

    # Checking that column names work.
    sce10x3 <- read10xCounts(tmpdir, col.names=TRUE)
    expect_identical(colnames(sce10x3), sce10x3$Barcode)

    sce10x4 <- read10xCounts(c(tmpdir, tmpdir), col.names=TRUE)
    expect_identical(colnames(sce10x4), paste0(rep(1:2, each=ncol(sce10x3)), "_", colnames(sce10x3)))

    # Checking that sample names work.
    sce10x5 <- read10xCounts(c(B=tmpdir, C=tmpdir))
    expect_identical(colData(sce10x5)$Sample, rep(c("B", "C"), each=ncol(sce10x3)))
    expect_identical(metadata(sce10x5)$Samples, c(B=tmpdir, C=tmpdir)) 

    sce10x6 <- read10xCounts(c(tmpdir, tmpdir), sample.names=c("A", "B"))
    expect_identical(colData(sce10x6)$Sample, rep(c("A", "B"), each=ncol(sce10x3)))
    expect_identical(metadata(sce10x6)$Samples, c(tmpdir, tmpdir)) 
})

test_that("read10xCounts works for sparse counts with odd inputs", {
    tmpdir <- tempfile()
    gene.symb2 <- paste0(gene.symb, sample(c("#", "'", '"', ""), length(gene.ids), replace=TRUE)) # full of weird elements.
    write10xCounts(path=tmpdir, my.counts, gene.id=gene.ids, gene.symbol=gene.symb2, barcodes=cell.ids)
    sce10x <- read10xCounts(tmpdir)

    expect_identical(assay(sce10x, withDimnames=FALSE), my.counts)
    expect_identical(colData(sce10x)$Barcode, cell.ids)
    expect_identical(rowData(sce10x)$ID, gene.ids)
    expect_identical(rowData(sce10x)$Symbol, gene.symb2)
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

    # Works in delayed mode.
    sce10delayed <- read10xCounts(c(tmpdir, tmpdir), delayed=TRUE)
    expect_s4_class(counts(sce10delayed), "DelayedMatrix")
    converted <- as(counts(sce10delayed), "dgCMatrix")
    expect_identical(converted, cbind(alt.counts, alt.counts))
})

test_that("read10xCounts works correctly for zipped files", {
    # Works for version 2:
    tmpdir <- tempfile()
    write10xCounts(path=tmpdir, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids)

    ref <- read10xCounts(tmpdir)
    lapply(list.files(tmpdir, full.names=TRUE), R.utils::gzip)
    alt <- read10xCounts(tmpdir)
    expect_identical(ref, alt)

    # Works for version 3:
    tmpdir <- tempfile()
    write10xCounts(path=tmpdir, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids, version="3")

    ref <- read10xCounts(tmpdir)
    lapply(list.files(tmpdir, full.names=TRUE), R.utils::gunzip)
    alt <- read10xCounts(tmpdir)
    expect_identical(ref, alt)
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

test_that("read10xCounts works correctly for prefixes", {
    tmpdir1 <- tempfile()
    write10xCounts(path=tmpdir1, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids)
    tmpdir2 <- tempfile()
    write10xCounts(path=tmpdir2, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids)

    tmpdir.all <- tempfile()
    dir.create(tmpdir.all, showWarnings=FALSE)
    first.files <- list.files(tmpdir1, full.names=TRUE)
    file.copy(first.files, file.path(tmpdir.all, paste0("jelly_", basename(first.files))))
    second.files <- list.files(tmpdir2, full.names=TRUE)
    file.copy(second.files, file.path(tmpdir.all, paste0("peanut_", basename(second.files))))

    out <- read10xCounts(file.path(tmpdir.all, c("jelly_", "peanut_")), type="prefix")
    alt <- read10xCounts(c(tmpdir1, tmpdir2))
    expect_identical(assay(out), assay(alt))
    expect_identical(out$Barcode, alt$Barcode)
})

test_that("read10xCounts works correctly with mismatching features", {
    tmpdir1 <- tempfile()
    write10xCounts(path=tmpdir1, my.counts, gene.id=gene.ids, gene.symbol=gene.symb, barcodes=cell.ids)

    tmpdir2 <- tempfile()
    keep <- 5:19
    write10xCounts(path=tmpdir2, my.counts[keep,], gene.id=gene.ids[keep], gene.symbol=gene.symb[keep], barcodes=cell.ids)

    expect_error(read10xCounts(c(tmpdir1, tmpdir2)), "gene information differs")

    # Intersection works as expected.
    sce10x <- read10xCounts(c(tmpdir1, tmpdir2), intersect.genes=TRUE)
    expect_identical(rownames(sce10x), gene.ids[keep])
    expect_identical(assay(sce10x, withDimnames=FALSE), cbind(my.counts[keep,], my.counts[keep,]))
})
