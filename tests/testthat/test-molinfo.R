# This tests that swappedDrops works correctly.
# library(DropletUtils); library(testthat); source("test-molinfo.R")

tmpdir <- tempfile()
dir.create(tmpdir)
ngenes <- 20L
barcode <- 4L

set.seed(909)
test_that("Extraction of molecule information fields works correctly", {
    output <- sim10xMolInfo(tmpdir, return.tab=TRUE, barcode=barcode, nsamples=3)
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

set.seed(908)
test_that("Automatic detection of the molecule information fields works correctly", {
    # Checking automatic detection of the barcode length works.
    for (blen in c(4L, 6L, 8L)) { 
        output <- sim10xMolInfo(tmpdir, barcode=blen)
        current <- read10xMolInfo(output)
        expect_true(all(nchar(current$data$cell)==blen))
    }
})


