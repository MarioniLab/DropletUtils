# This tests that extraction of molecule information works correctly.
# library(DropletUtils); library(testthat); source("test-molinfo.R")

tmpdir <- tempfile(fileext=".h5")
barcode <- 4L
ngenes <- 20L

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
    output <- DropletUtils:::simBasicMolInfo(tmpdir, return.tab=TRUE, barcode=barcode)

    combined <- output$original
    combined <- combined[combined$gene<=ngenes,]

    current <- read10xMolInfo(output$files, barcode.length=barcode)
    expect_identical(length(current$genes), ngenes)
    expect_identical(as.integer(current$data$umi), combined$umi)
    expect_identical(as.integer(current$data$gene), combined$gene)
    expect_identical(as.integer(current$data$reads), combined$reads)
    expect_identical(as.integer(current$data$gem_group), rep(1L, nrow(combined)))

    # Checking that there is a 1:1 relationship between the cell barcodes and cell IDs.
    by.barcode <- split(combined$cell, current$data$cell)
    expect_true(all(lengths(lapply(by.barcode, unique))==1L))
    by.cell.id <- split(current$data$cell, combined$cell)
    expect_true(all(lengths(lapply(by.cell.id, unique))==1L))

    # Checking that using too little barcode length underestimates the number of cells.
    current2 <- read10xMolInfo(output$files, barcode.length=barcode - 1L)
    expect_true(length(unique(current2$data$cell)) < length(unique(current$data$cell)))
    current3 <- read10xMolInfo(output$files, barcode.length=barcode + 1L)
    expect_true(length(unique(current3$data$cell)) == length(unique(current$data$cell)))
})

set.seed(9091)
test_that("Extraction of subsets of the molinfo fields works correctly", {
    output <- DropletUtils:::simBasicMolInfo(tmpdir, barcode=barcode)
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

    subbed <- read10xMolInfo(output, get.gem=FALSE, get.reads=FALSE, get.cell=FALSE, get.gene=FALSE, get.umi=FALSE, 
        get.library=FALSE, keep.unmapped=TRUE)
    expect_identical(nrow(subbed$data), nrow(fullun$data))
    expect_identical(ncol(subbed$data), 0L)
    expect_identical(subbed$genes, fullun$genes)
})

set.seed(908)
test_that("Automatic detection of the molecule information fields works correctly", {
    for (blen in c(4L, 6L, 8L)) { 
        output <- DropletUtils:::simBasicMolInfo(tmpdir, barcode=blen)
        current <- read10xMolInfo(output)
        expect_true(all(nchar(current$data$cell)==blen))
    }
})

set.seed(908)
test_that("read10xMolInfo responds correctly to the CellRanger version", {
    set.seed(100)
    output <- DropletUtils:::simBasicMolInfo(tmpdir, barcode=6, version="2")
    expect_false("barcode_idx" %in% rhdf5::h5ls(output)$name)
    restored <-read10xMolInfo(output)

    set.seed(100)
    output <- DropletUtils:::simBasicMolInfo(tmpdir, barcode=6, version="3")
    expect_true("barcode_idx" %in% rhdf5::h5ls(output)$name)
    restored2 <-read10xMolInfo(output)

    expect_type(restored2$data$library, "integer")
    restored2$data$library <- NULL
    expect_identical(restored$data, restored2$data)

    # Pulling out the library information as JSON.
    full.plus <- read10xMolInfo(output, extract.library.info=TRUE)
    expect_identical(length(full.plus$library.info), 3L)
})

set.seed(907)
test_that("read10xMolInfo works with silly inputs containing no molecules", {
    out.paths <- DropletUtils:::simBasicMolInfo(tmpdir, nmolecules=0, ngenes=ngenes, barcode=barcode)
    out <- read10xMolInfo(out.paths, barcode=barcode)
    expect_identical(nrow(out$data), 0L)
    expect_identical(length(out$genes), ngenes)

    # Checking that it doesn't throw up with automatic barcode detection.
    out2 <- read10xMolInfo(out.paths)
    expect_identical(out, out2)
   
    # Checking that it behaves when there aren't even any genes. 
    out.paths <- DropletUtils:::simBasicMolInfo(tmpdir, nmolecules=0, ngenes=0, barcode.length=barcode) 
    out <- read10xMolInfo(out.paths, barcode=barcode)
    expect_identical(nrow(out$data), 0L)
    expect_identical(length(out$genes), 0L)
})

set.seed(908)
test_that("Molecule information extraction by library works correctly", {
    out.path <- DropletUtils:::simBasicMolInfo(tmpdir, barcode=barcode, version="3")

    # Naive extraction.
    ref <- read10xMolInfo(out.path, get.library=FALSE)
    out <- DropletUtils:::.extract_mol_info(out.path)
    expect_identical(out, ref)

    # Subsetting to the first library. 
    ref <- read10xMolInfo(out.path, get.library=TRUE)

    out <- DropletUtils:::.extract_mol_info(out.path, use.library=1)
    sub <- ref$data[ref$data$library == 1,]
    sub$library <- NULL
    expect_identical(out$data, sub)

    out <- DropletUtils:::.extract_mol_info(out.path, use.library=c(3, 1))
    sub <- ref$data[ref$data$library %in% c(1, 3),]
    sub$library <- NULL
    expect_identical(out$data, sub)

    # Works with character vectors.
    out <- DropletUtils:::.extract_mol_info(out.path, use.library="B")
    sub <- ref$data[ref$data$library == 2,]
    sub$library <- NULL
    expect_identical(out$data, sub)

    out <- DropletUtils:::.extract_mol_info(out.path, use.library=c("A", "B"))
    sub <- ref$data[ref$data$library %in% 1:2,]
    sub$library <- NULL
    expect_identical(out$data, sub)

    # Works with subsetting.
    out <- DropletUtils:::.extract_mol_info(out.path, use.library=2, subset.library.features=TRUE)
    sub <- ref$data[ref$data$library == 2,]
    sub$gene <- match(sub$gene, which(ref$feature.type=="B"))
    sub$library <- NULL
    expect_identical(out$data, sub)
    expect_identical(out$genes, ref$genes[ref$feature.type=="B"])

    out <- DropletUtils:::.extract_mol_info(out.path, use.library="C", subset.library.features=TRUE)
    sub <- ref$data[ref$data$library == 3,]
    sub$gene <- match(sub$gene, which(ref$feature.type=="C"))
    sub$library <- NULL
    expect_identical(out$data, sub)
    expect_identical(out$genes, ref$genes[ref$feature.type=="C"])

    # Preserves library information if requested.
    out <- DropletUtils:::.extract_mol_info(out.path, use.library="A", get.library=TRUE, extract.library.info=TRUE)
    ref2 <- read10xMolInfo(out.path, extract.library.info=TRUE)
    ref2$data <- ref2$data[ref2$data$library == 1,]
    expect_identical(out, ref2)
})
