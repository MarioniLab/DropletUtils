# This tests that extraction of molecule information works correctly.
# library(DropletUtils); library(testthat); source("test-molinfo.R")

tmpdir <- tempfile()
dir.create(tmpdir)
ngenes <- 20L
barcode <- 4L

set.seed(910)
test_that("barcode extraction is working correctly", {
    library(rhdf5)
    for (blen in c(4, 7, 10)) {
         all.barcodes <- sample(4^blen, 10000, replace=TRUE) - 1L

         out.file <- tempfile(fileext="h5")
         h5 <- h5createFile(out.file)
         h5write(all.barcodes, out.file, "barcode")

         out <- DropletUtils:::get_cell_barcodes(out.file, "barcode", blen)
         guess <- DropletUtils:::get_cell_barcodes(out.file, "barcode", NULL)
         expect_identical(out, guess)

         # Manually doing the bit masks.
         progressive <- DropletUtils:::.unmask_barcode(all.barcodes, blen)
         expect_identical(out, progressive)
    }
})

set.seed(909)
test_that("Extraction of molecule information fields works correctly", {
    output <- DropletUtils:::sim10xMolInfo(tmpdir, return.tab=TRUE, barcode=barcode, nsamples=3)
    for (i in seq_along(output$files)) {
        ref.original <- output$original[output$original$sample==i,]        
        ref.swapped <- output$swapped[output$swapped$sample==i,]
        combined <- rbind(ref.original, ref.swapped)
        combined <- combined[combined$gene<ngenes,]

        current <- read10xMolInfo(output$files[i], barcode.length=barcode)
        expect_identical(length(current$genes), ngenes)
        expect_identical(as.integer(current$data$umi), combined$umi)
        expect_identical(as.integer(current$data$gene), combined$gene+1L)
        expect_identical(as.integer(current$data$reads), combined$reads)
        expect_identical(as.integer(current$data$gem_group), rep(1L, nrow(combined)))

        # Checking that there is a 1:1 relationship between the cell barcodes and cell IDs.
        by.barcode <- split(combined$cell, current$data$cell)
        expect_true(all(lengths(lapply(by.barcode, unique))==1L))
        by.cell.id <- split(current$data$cell, combined$cell)
        expect_true(all(lengths(lapply(by.cell.id, unique))==1L))

        # Checking that using too little barcode length underestimates the number of cells.
        current2 <- read10xMolInfo(output$files[i], barcode.length=barcode - 1L)
        expect_true(length(unique(current2$data$cell)) < length(unique(current$data$cell)))
        current3 <- read10xMolInfo(output$files[i], barcode.length=barcode + 1L)
        expect_true(length(unique(current3$data$cell)) == length(unique(current$data$cell)))
    }
})

set.seed(9091)
test_that("Extraction of subsets of the molinfo fields works correctly", {
    output <- DropletUtils:::sim10xMolInfo(tmpdir, barcode=barcode, nsamples=1)
    full <- read10xMolInfo(output)
    
    # Discounting the GEM.
    subbed <- read10xMolInfo(output, get.gem=FALSE)
    tmp <- full
    tmp$data$gem_group <- NULL
    expect_identical(subbed, tmp)

    # Discounting the genes. 
    subbed <- read10xMolInfo(output, get.gene=FALSE)
    tmp <- full
    tmp$data$gene <- NULL
    expect_identical(subbed, tmp)

    fullun <- read10xMolInfo(output, keep.unmapped=TRUE)
    subbed <- read10xMolInfo(output, get.gene=FALSE, keep.unmapped=TRUE)
    fullun$data$gene <- NULL
    expect_identical(subbed, fullun)

    # Discounting everything.
    subbed <- read10xMolInfo(output, get.gem=FALSE, get.reads=FALSE, get.cell=FALSE, get.gene=FALSE, get.umi=FALSE)
    expect_identical(nrow(subbed$data), nrow(full$data))
    expect_identical(ncol(subbed$data), 0L)
    expect_identical(subbed$genes, full$genes)

    subbed <- read10xMolInfo(output, get.gem=FALSE, get.reads=FALSE, get.cell=FALSE, get.gene=FALSE, get.umi=FALSE, keep.unmapped=TRUE)
    expect_identical(nrow(subbed$data), nrow(fullun$data))
    expect_identical(ncol(subbed$data), 0L)
    expect_identical(subbed$genes, fullun$genes)
})

set.seed(908)
test_that("Automatic detection of the molecule information fields works correctly", {
    for (blen in c(4L, 6L, 8L)) { 
        output <- DropletUtils:::sim10xMolInfo(tmpdir, barcode=blen)
        current <- read10xMolInfo(output)
        expect_true(all(nchar(current$data$cell)==blen))
    }
})

set.seed(908)
test_that("read10xMolInfo responds correctly to the CellRanger version", {
    set.seed(100)
    output <- DropletUtils:::sim10xMolInfo(tmpdir, barcode=6, version="2")
    expect_false("barcode_idx" %in% rhdf5::h5ls(output)$name)
    restored <-read10xMolInfo(output)

    set.seed(100)
    output <- DropletUtils:::sim10xMolInfo(tmpdir, barcode=6, version="3")
    expect_true("barcode_idx" %in% rhdf5::h5ls(output)$name)
    restored2 <-read10xMolInfo(output)

    expect_identical(restored, restored2)
})

set.seed(907)
test_that("read10xMolInfo works with silly inputs containing no molecules", {
    out.paths <- DropletUtils:::sim10xMolInfo(tmpdir, nsamples=1, nmolecules=0, swap.frac=0, ngenes=ngenes, barcode=barcode)
    out <- read10xMolInfo(out.paths, barcode=barcode)
    expect_identical(nrow(out$data), 0L)
    expect_identical(length(out$genes), ngenes)

    # Checking that it doesn't throw up with automatic barcode detection.
    out2 <- read10xMolInfo(out.paths)
    expect_identical(out, out2)
   
    # Checking  that it behaves when there aren't even any genes. 
    out.paths <- DropletUtils:::sim10xMolInfo(tmpdir, nsamples=1, nmolecules=0, ngenes=0, swap.frac=0, barcode.length=barcode) 
    out <- read10xMolInfo(out.paths, barcode=barcode)
    expect_identical(nrow(out$data), 0L)
    expect_identical(length(out$genes), 0L)
})

