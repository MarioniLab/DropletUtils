# This tests the ability to ingest various 10x matrix files.
# library(testthat); library(DropletUtils); source("test-read10x.R")

library(DropletTestFiles)

test_that("read10xCounts works for version 2 matrices (tarball)", {
    # Filtered:
    fname <- getTestFile("tenx-2.1.0-pbmc4k/1.0.0/filtered.tar.gz")
    tmp <- tempfile()
    untar(fname, exdir=tmp)
    out <- read10xCounts(file.path(tmp, "filtered_gene_bc_matrices/GRCh38"))

    expect_false(is.null(rownames(out)))
    expect_type(out$Barcode, "character")
    expect_s4_class(counts(out), "CsparseMatrix")

    # Raw:
    fname <- getTestFile("tenx-2.1.0-pbmc4k/1.0.0/raw.tar.gz")
    tmp <- tempfile()
    untar(fname, exdir=tmp)
    out2 <- read10xCounts(file.path(tmp, "raw_gene_bc_matrices/GRCh38"))

    expect_false(is.null(rownames(out2)))
    expect_type(out2$Barcode, "character")
    expect_s4_class(counts(out2), "CsparseMatrix")

    # Comparing the two.
    expect_identical(rowData(out), rowData(out2))
    expect_true(ncol(out2) > ncol(out))
    expect_true(all(out$Barcode %in% out2$Barcode))
})

test_that("read10xCounts works for version 2 matrices (HDF5)", {
    fname <- getTestFile("tenx-2.1.0-pbmc4k/1.0.0/raw.h5")
    out <- read10xCounts(fname, type="HDF5")

    expect_false(is.null(rownames(out)))
    expect_type(out$Barcode, "character")
    expect_s4_class(counts(out), "DelayedMatrix")
    expect_type(counts(out)[1,], "integer")
})

test_that("read10xCounts works for version 3 matrices (tarball)", {
    # Filtered:
    fname <- getTestFile("tenx-3.1.0-5k_pbmc_protein_v3/1.0.0/filtered.tar.gz")
    tmp <- tempfile()
    untar(fname, exdir=tmp)
    out <- read10xCounts(file.path(tmp, "filtered_feature_bc_matrix"))

    expect_identical(length(unique(rowData(out)$Type)), 2L)
    expect_false(is.null(rownames(out)))
    expect_type(out$Barcode, "character")
    expect_s4_class(counts(out), "CsparseMatrix")

    # Raw:
    fname <- getTestFile("tenx-3.1.0-5k_pbmc_protein_v3/1.0.0/raw.tar.gz")
    tmp <- tempfile()
    untar(fname, exdir=tmp)
    out2 <- read10xCounts(file.path(tmp, "raw_feature_bc_matrix"))

    expect_identical(length(unique(rowData(out2)$Type)), 2L)
    expect_false(is.null(rownames(out2)))
    expect_type(out2$Barcode, "character")
    expect_s4_class(counts(out2), "CsparseMatrix")

    # Comparing the two.
    expect_identical(rowData(out), rowData(out2))
    expect_true(ncol(out2) > ncol(out))
    expect_true(all(out$Barcode %in% out2$Barcode))
})

test_that("read10xCounts works for version 3 matrices (HDF5)", {
    # Filtered:
    fname <- getTestFile("tenx-3.1.0-5k_pbmc_protein_v3/1.0.0/filtered.h5")
    out <- read10xCounts(fname, type="HDF5")

    expect_identical(length(unique(rowData(out)$Type)), 2L)
    expect_false(is.null(rownames(out)))
    expect_type(out$Barcode, "character")
    expect_s4_class(counts(out), "DelayedMatrix")
    expect_type(counts(out)[1,], "integer")

    # Raw:
    fname <- getTestFile("tenx-3.1.0-5k_pbmc_protein_v3/1.0.0/raw.h5")
    out2 <- read10xCounts(fname, type="HDF5")

    expect_identical(length(unique(rowData(out2)$Type)), 2L)
    expect_false(is.null(rownames(out2)))
    expect_type(out2$Barcode, "character")
    expect_s4_class(counts(out2), "DelayedMatrix")
    expect_type(counts(out2)[1,], "integer")

    # Comparing the two.
    expect_identical(rowData(out), rowData(out2))
    expect_true(ncol(out2) > ncol(out))
    expect_true(all(out$Barcode %in% out2$Barcode))
})

test_that("read10xCounts works for version 4 matrices (tarball)", {
    # Filtered:
    fname <- getTestFile("tenx-4.0.0-SC3_v3_NextGem_DI_Neuron_10K/1.0.0/filtered.tar.gz")
    tmp <- tempfile()
    untar(fname, exdir=tmp)

    out <- read10xCounts(file.path(tmp, "filtered_feature_bc_matrix"))

    expect_identical(length(unique(rowData(out)$Type)), 1L)
    expect_false(is.null(rownames(out)))
    expect_type(out$Barcode, "character")
    expect_s4_class(counts(out), "CsparseMatrix")

    # Raw:
    fname <- getTestFile("tenx-4.0.0-SC3_v3_NextGem_DI_Neuron_10K/1.0.0/raw.tar.gz")
    tmp <- tempfile()
    untar(fname, exdir=tmp)

    out2 <- read10xCounts(file.path(tmp, "raw_feature_bc_matrix"))

    expect_identical(length(unique(rowData(out)$Type)), 1L)
    expect_false(is.null(rownames(out)))
    expect_type(out$Barcode, "character")
    expect_s4_class(counts(out), "CsparseMatrix")

    # Comparing the two.
    expect_identical(rowData(out), rowData(out2))
    expect_true(ncol(out2) > ncol(out))
    expect_true(all(out$Barcode %in% out2$Barcode))
})

test_that("read10xCounts works for version 4 matrices (HDF5)", {
    # Filtered:
    fname <- getTestFile("tenx-4.0.0-SC3_v3_NextGem_DI_Neuron_10K/1.0.0/filtered.h5")
    out <- read10xCounts(fname, type="HDF5")

    expect_identical(length(unique(rowData(out)$Type)), 1L)
    expect_false(is.null(rownames(out)))
    expect_type(out$Barcode, "character")
    expect_s4_class(counts(out), "DelayedMatrix")
    expect_type(counts(out)[1,], "integer")

    # Raw:
    fname <- getTestFile("tenx-4.0.0-SC3_v3_NextGem_DI_Neuron_10K/1.0.0/raw.h5")
    out2 <- read10xCounts(fname, type="HDF5")

    expect_identical(length(unique(rowData(out2)$Type)), 1L)
    expect_false(is.null(rownames(out2)))
    expect_type(out2$Barcode, "character")
    expect_s4_class(counts(out2), "DelayedMatrix")
    expect_type(counts(out2)[1,], "integer")

    # Comparing the two.
    expect_identical(rowData(out), rowData(out2))
    expect_true(ncol(out2) > ncol(out))
    expect_true(all(out$Barcode %in% out2$Barcode))
})
