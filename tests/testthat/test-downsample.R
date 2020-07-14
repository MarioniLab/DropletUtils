# Testing the downsampleReads function.
# library(DropletUtils); library(testthat); source("test-downsample.R")

library(Matrix)
set.seed(501)

test_that("downsampling from the reads yields correct results", {
    barcode <- 4L
    tmpdir <- tempfile()
    dir.create(tmpdir)
    out.paths <- DropletUtils:::simBasicMolInfo(tmpdir, ngenes=100, barcode.length=barcode) 
   
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

    # Making it easier to check the totals, by making all UMIs have a read count of 1.
    out.paths <- DropletUtils:::simBasicMolInfo(tmpdir, ngenes=100, barcode.len=barcode, ave.read=0) 
    full.tab <- downsampleReads(out.paths, barcode, prop=1)
    expect_equal(sum(downsampleReads(out.paths, barcode, prop=0.555)), round(0.555*sum(full.tab))) # Again, avoiding rounding differences.
    expect_equal(sum(downsampleReads(out.paths, barcode, prop=0.111)), round(0.111*sum(full.tab)))
    expect_equal(colSums(downsampleReads(out.paths, barcode, prop=0.555, bycol=TRUE)), round(0.555*colSums(full.tab)))
    expect_equal(colSums(downsampleReads(out.paths, barcode, prop=0.111, bycol=TRUE)), round(0.111*colSums(full.tab)))

    # Checking behaviour on silly inputs where there are no reads, or no genes.
    ngenes <- 20L
    out.paths <- DropletUtils:::simBasicMolInfo(tmpdir, nmolecules=0, ngenes=ngenes, barcode.length=barcode) 
    out <- downsampleReads(out.paths, barcode, prop=0.5)
    expect_identical(dim(out), c(ngenes, 0L))

    out.paths <- DropletUtils:::simBasicMolInfo(tmpdir, nmolecules=0, ngenes=0, barcode.length=barcode) 
    out <- downsampleReads(out.paths, barcode, prop=0.5)
    expect_identical(dim(out), c(0L, 0L))
})
    
test_that("downsampling from the reads compares correctly to downsampleMatrix", {
    # Manually creating files for comparison to downsampleMatrix - this relies on ordered 'gene' and 'cell', 
    # so that the retention probabilities applied to each molecule are the same across functions.
    ngenes <- 4
    gene.count <- seq_len(ngenes)*100
    ncells <- 200
    nmolecules <- sum(gene.count)*ncells

    tmpdir <- tempfile()
    dir.create(tmpdir)
    out.file <- file.path(tmpdir, "out.h5")

    library(rhdf5)
    h5 <- h5createFile(out.file)
    h5write(rep(seq_len(ncells), each=sum(gene.count)), out.file, "barcode")
    h5write(seq_len(nmolecules), out.file, "umi")
    h5write(rep(rep(seq_len(ngenes)-1L, gene.count), ncells), out.file, "gene")
    h5write(rep(1, nmolecules), out.file, "gem_group")
    h5write(rep(1, nmolecules), out.file, "reads") # one read per molecule.
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
})

test_that("downsampling from the reads works correctly with feature subsets", {
    barcode <- 4L
    tmpdir <- tempfile()
    dir.create(tmpdir)
    out.paths <- DropletUtils:::simBasicMolInfo(tmpdir, ngenes=100, barcode.length=barcode) 

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
