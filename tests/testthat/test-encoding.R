# Tests the encodeSequences function works correctly.
# library(DropletUtils); library(testthat); source("test-encoding.R")

CONVERTER <- function(sequences) {
    as.seq <- strsplit(sequences, "")
    charmat <-  do.call(rbind, as.seq)
    idx.out <- integer(length(as.seq))
    mult <-  1L
    for (i in ncol(charmat):1) {
        idx.out <- idx.out + (match(charmat[,i], c("A", "C", "G", "T")) - 1L) * mult
        mult <- mult * 4L
    } 
    return(idx.out)
}

test_that("encodeSequences works as expected", {
    stuff <- c("ACACGTGT", "GCGATCAA", "TTCAACGG", "GGGGGGGG")
    expect_identical(encodeSequences(stuff), CONVERTER(stuff))
    expect_identical(encodeSequences(strrep("A", 1:10)), integer(10))
})

test_that("encodeSequences behaves with silly inputs", {
    expect_identical(encodeSequences(character(0)), integer(0))
    expect_error(encodeSequences("banana"), "banana")
    expect_error(encodeSequences("ACAGCTAGCATCGATC"), "15 nt")
})
