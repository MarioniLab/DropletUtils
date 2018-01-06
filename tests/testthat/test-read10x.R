# Testing the read10xCounts function.
# library(DropletUtils); library(testthat); source("test-read10x.R")

test_that("read10xCounts works correctly", {
    # Mocking up some 10X genomics output.
    tmpdir <- tempfile()
    dir.create(tmpdir)

    library(Matrix)
    my.counts <- matrix(rpois(1000, lambda=5), ncol=10, nrow=100)
    my.counts <- as(my.counts, "dgCMatrix")
    writeMM(my.counts, file=file.path(tmpdir, "matrix.mtx"))

    ngenes <- nrow(my.counts)
    gene.ids <- paste0("GENE", seq_len(ngenes))
    gene.symb <- paste0(sample(LETTERS, replace=TRUE, ngenes),
                        sample(LETTERS, replace=TRUE, ngenes),
                        sample(LETTERS, replace=TRUE, ngenes), "-",
                        sample(9, replace=TRUE, ngenes))
    write.table(data.frame(gene.ids, gene.symb), file=file.path(tmpdir, "genes.tsv"),
                    quote=FALSE, col.names=FALSE, row.names=FALSE)   

    cell.ids <- paste0("BARCODE-", seq_len(ncol(my.counts)))
    write(cell.ids, file=file.path(tmpdir, "barcodes.tsv"))

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
