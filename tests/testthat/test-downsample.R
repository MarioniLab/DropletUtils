# Testing the downsampleReads function.
# library(DropletUtils); library(testthat); source("test-downsample.R")

library(Matrix)
set.seed(501)
barcode <- 4L

test_that("downsampling from the reads yields correct results", {
    tmpfile <- tempfile(fileext=".h5")
    out.paths <- DropletUtils:::simBasicMolInfo(tmpfile, ngenes=100, barcode.length=barcode) 
   
    # Creating the full matrix, and checking that it's the same when no downsampling is requested.
    collated <- read10xMolInfo(out.paths, barcode)
    all.cells <- sort(unique(collated$data$cell))
    full.tab <- makeCountMatrix(collated$data$gene, collated$data$cell, all.genes=collated$genes)
    colnames(full.tab) <- paste0(colnames(full.tab), "-1")

    out <- downsampleReads(out.paths, barcode, prop=1)
    expect_equal(out, full.tab)

    # Checking that the ordering of cells is equivalent.
    stats <- get10xMolInfoStats(out.paths)
    expect_identical(colnames(out), sprintf("%s-%i", stats$cell, stats$gem_group))

    # Checking that some downsampling has occurred (hard to check the totals, as UMI counts != read counts).
    for (down in 1:4/11) {
        out <- downsampleReads(out.paths, barcode, prop=down)
        expect_true(all(out <= full.tab))
        expect_false(all(out==full.tab))
    }
})

test_that("downsampling from the reads yields correct results", {
    tmpfile <- tempfile(fileext=".h5")

    # Making it easier to check the totals, by making all UMIs have a read count of 1.
    out.paths <- DropletUtils:::simBasicMolInfo(tmpfile, ngenes=100, barcode.len=barcode, ave.read=0) 
    full.tab <- downsampleReads(out.paths, barcode, prop=1)

    # Setting an odd fraction to avoid rounding differences between OS's.
    expect_equal(sum(downsampleReads(out.paths, barcode, prop=0.555)), round(0.555*sum(full.tab))) 
    expect_equal(sum(downsampleReads(out.paths, barcode, prop=0.111)), round(0.111*sum(full.tab)))
    expect_equal(colSums(downsampleReads(out.paths, barcode, prop=0.555, bycol=TRUE)), round(0.555*colSums(full.tab)))
    expect_equal(colSums(downsampleReads(out.paths, barcode, prop=0.111, bycol=TRUE)), round(0.111*colSums(full.tab)))

    # Checking behaviour on silly inputs where there are no reads, or no genes.
    ngenes <- 20L
    out.paths <- DropletUtils:::simBasicMolInfo(tmpfile, nmolecules=0, ngenes=ngenes, barcode.length=barcode) 
    out <- downsampleReads(out.paths, barcode, prop=0.5)
    expect_identical(dim(out), c(ngenes, 0L))

    out.paths <- DropletUtils:::simBasicMolInfo(tmpfile, nmolecules=0, ngenes=0, barcode.length=barcode) 
    out <- downsampleReads(out.paths, barcode, prop=0.5)
    expect_identical(dim(out), c(0L, 0L))
})

test_that("downsampling from the reads works correctly with feature subsets", {
    tmpfile <- tempfile(fileext=".h5")
    out.paths <- DropletUtils:::simBasicMolInfo(tmpfile, ngenes=100, barcode.length=barcode) 

    # Full matrix is correctly extracted without any downsampling.
    collated <- read10xMolInfo(out.paths, barcode)
    all.cells <- sort(unique(collated$data$cell))
    full.tab <- makeCountMatrix(collated$data$gene, collated$data$cell, all.genes=collated$genes)
    colnames(full.tab) <- paste0(colnames(full.tab), "-1")

    features <- collated$genes[sort(sample(length(collated$genes), length(collated$genes)/2))]
    out <- downsampleReads(out.paths, barcode, prop=1, features=features)
    expect_equal(out, full.tab[as.character(features),])
    expect_identical(rownames(out), features)

    # Downsampling behaves as expected.
    out2 <- downsampleReads(out.paths, barcode, prop=0.1, features=features)
    expect_identical(rownames(out), rownames(out2))
    expect_true(all(out >= out2))
    expect_true(sum(out) >= sum(out2))
})

test_that("downsampling from the reads works correctly with library subsets", {
    tmpfile <- tempfile(fileext=".h5")
    out.paths <- DropletUtils:::simBasicMolInfo(tmpfile, ngenes=100, barcode.length=barcode, version="3")

    set.seed(100)
    out2 <- downsampleReads(out.paths, barcode, prop=0.1, use.library="A")

    info <- read10xMolInfo(tmpfile)
    expect_identical(rownames(out2), info$genes[info$feature.type=="A"])

    set.seed(100)
    ref <- downsampleReads(out.paths, barcode, prop=0.1, features=rownames(out2))
    expect_identical(out2, ref)
})
