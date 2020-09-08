# We can read in molecule information correctly from real datasets.
# library(DropletUtils); library(testthat); source("test-molinfo.R")

library(DropletTestFiles)

test_that("works for version 2", {
    fname <- getTestFile("tenx-2.1.0-pbmc4k/1.0.0/mol_info.h5")
    out <- read10xMolInfo(fname, extract.library.info=TRUE)

    expect_s4_class(out$data, "DataFrame")
    expect_true(nrow(out$data) > 1e7)
    expect_true(all(out$data$gene <= length(out$genes)+1L))

    # Clearing it out for the next element.
    rm(out)
    gc()
})

test_that("works for version 3", {
    fname <- getTestFile("tenx-3.0.0-pbmc_10k_protein_v3/1.0.0/mol_info.h5")
    out <- read10xMolInfo(fname, extract.library.info=TRUE)

    expect_s4_class(out$data, "DataFrame")
    expect_true(nrow(out$data) > 1e8)

    expect_true(all(out$data$gene <= length(out$genes)+1L))
    expect_identical(length(out$genes), length(out$feature.type))
    expect_identical(length(unique(out$feature.type)), 2L)

    expect_identical(length(out$library.info), 2L)
    expect_true(all(out$feature.type %in% vapply(out$library.info, "[[", name="library_type", "")))

    # Clearing it out for the next element.
    rm(out)
    gc()
})

test_that("works for version 4", {
    fname <- getTestFile("tenx-4.0.0-SC3_v3_NextGem_DI_Neuron_10K/1.0.0/mol_info.h5")
    out <- read10xMolInfo(fname, extract.library.info=TRUE)

    expect_s4_class(out$data, "DataFrame")
    expect_true(nrow(out$data) > 1e8)

    expect_true(all(out$data$gene <= length(out$genes)+1L))
    expect_identical(length(out$genes), length(out$feature.type))
    expect_identical(length(unique(out$feature.type)), 1L)

    expect_identical(length(out$library.info), 1L)
    expect_true(all(out$feature.type %in% vapply(out$library.info, "[[", name="library_type", "")))

    # Clearing it out for the next element.
    rm(out)
    gc()
})
