# This tests the chimericDrops functionality.
# library(DropletUtils); library(testthat); source("test-chimera.R")

set.seed(100001)

EDIT_H5 <- function(sample, subset) {
    umis <- rhdf5::h5read(sample, "umi")
    new.umis <- sample(head(umis, subset), length(umis), replace=TRUE) 
    rhdf5::h5write(new.umis, sample, "umi")
}

test_that("chimericDrops removes no duplicated UMIs with standard simulations", {
    tmpdir <- tempfile()

    # Default is to not remove anything,
    # as UMI sampling is done without replacement.
    fname <- DropletUtils:::simBasicMolInfo(tmpdir)

    output <- chimericDrops(fname, get.chimeric=TRUE)

    info <- read10xMolInfo(fname)
    ref <- makeCountMatrix(info$data$gene, info$data$cell, all.genes=info$genes)
    expect_identical(output$cleaned, ref)
    expect_identical(sum(abs(output$chimeric)), 0)    
})

test_that("chimericDrops removes duplicated UMIs in more complex cases", {
    for (actual.unique in c(1000, 5000, 9000)) {
        tmpdir <- tempfile()
        fname <- DropletUtils:::simBasicMolInfo(tmpdir)
        EDIT_H5(fname, subset=actual.unique) # introducing duplicates.

        info <- read10xMolInfo(fname)
        ref <- makeCountMatrix(info$data$gene, info$data$cell, all.genes=info$genes)

        output <- chimericDrops(fname, get.chimeric=TRUE)
        expect_identical(output$cleaned + output$chimeric, ref)

        # More sophisticated check, setting min.frac to basically remove all duplicates
        # and removing unmapped molecules for testing simplicity.
        keep <- info$data$gene!=length(info$genes)
        info2 <- info$data[keep,]

        output <- removeChimericDrops(
            cells=info2$cell, 
            umis=info2$umi,
            genes=info2$gene,
            nreads=info2$reads,
            ref.genes=info$genes,
            min.frac=1,
            get.chimeric=TRUE
        )

        ref <- makeCountMatrix(info2$gene, info2$cell, all.genes=info$genes)
        expect_identical(output$cleaned + output$chimeric, ref)

        by <- info2[,c("umi", "cell")]
        expect_equal(
            sum(table(selfmatch(by))==1),
            sum(output$cleaned)
        )
    }
})

test_that("chimericDrops respects the use.library= restriction", {
    tmpdir <- tempfile(fileext=".h5")
    fname <- DropletUtils:::simBasicMolInfo(tmpdir, version="3")

    # Behaves properly when no restriction is placed down.
    ref <- chimericDrops(fname, get.chimeric=TRUE)
    expect_true(all(dim(ref[[1]]) > 0L))
    expect_true(all(dim(ref[[2]]) > 0L))

    ref2 <- chimericDrops(fname, get.chimeric=TRUE, use.library=1:3)
    expect_identical(ref, ref2)

    # Correctly empties out when a restriction is applied.
    output <- chimericDrops(fname, get.chimeric=TRUE, use.library="XXX")
    expect_true(all(dim(output[[1]]) == 0L))
    expect_true(all(dim(output[[2]]) == 0L))

    output <- chimericDrops(fname, get.chimeric=TRUE, use.library="A")
    expect_true(nrow(output[[1]]) < nrow(ref[[1]]))
})

test_that("chimericDrops reports diagnostics correctly", {
    for (actual.unique in c(1000, 5000, 9000)) {
        tmpdir <- tempfile()
        fname <- DropletUtils:::simBasicMolInfo(tmpdir)
        EDIT_H5(fname, subset=actual.unique) # introducing duplicates.

        info <- read10xMolInfo(fname)
        output <- chimericDrops(fname, get.chimeric=TRUE, get.diagnostics=TRUE)
        by <- info$data[,c("umi", "cell")]
        expect_identical(sum(!duplicated(by)), nrow(output$diagnostics))
    }
})
