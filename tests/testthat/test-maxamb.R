# This tests the maximum ambience function.
# library(testthat); library(DropletUtils); source("setup.R"); source("test-maxamb.R")

# Making up some data.
set.seed(1100000)
ambient <- c(runif(900, 0, 0.1), runif(100))
y <- rpois(1000, ambient * 50)
y2 <- y + rpois(1000, 5) # adding some actual biology.

test_that("maximumAmbience function works as expected", {
    scaling <- maximumAmbience(y, ambient)
    expect_equal(scaling, sum(y)/sum(ambient))
    
    scaling <- maximumAmbience(y2, ambient)
    expect_true(scaling > 50)
    expect_true(scaling < sum(y2)/sum(ambient))
})

test_that("maximumAmbience function handles matrix ambience", {
    stuff <- maximumAmbience(cbind(y, y2), ambient)
    expect_identical(unname(stuff[1]), maximumAmbience(y, ambient))
    expect_identical(unname(stuff[2]), maximumAmbience(y2, ambient))

    ambient2 <- c(runif(900, 0, 0.1), runif(100))
    stuff <- maximumAmbience(cbind(y, y2), cbind(ambient, ambient2))
    expect_identical(unname(stuff[1]), maximumAmbience(y, ambient))
    expect_identical(unname(stuff[2]), maximumAmbience(y2, ambient2))

    stuff <- maximumAmbience(cbind(y, y2), cbind(ambient, ambient2), mode="profile")
    expect_identical(unname(stuff[,1,drop=FALSE]), unname(maximumAmbience(y, ambient, mode="profile")))
    expect_identical(unname(stuff[,2,drop=FALSE]), unname(maximumAmbience(y2, ambient2, mode="profile")))
})

test_that("maximumAmbience works correctly in parallel", {
    combined <- cbind(y, y2)
    ref <- maximumAmbience(y, ambient)
    scaling <- maximumAmbience(y, ambient, BPPARAM=safeBPPARAM)
    expect_identical(ref, scaling)

    ambient2 <- c(runif(900, 0, 0.1), runif(100))
    amb.mat <- cbind(ambient, ambient2)
    ref <- maximumAmbience(combined, amb.mat)
    scaling <- maximumAmbience(combined, amb.mat, BPPARAM=safeBPPARAM)
    expect_identical(ref, scaling)
})

test_that("maximumAmbience function handles dispersion and threshold changes", {
    scaling <- maximumAmbience(y2, ambient)

    scaling2 <- maximumAmbience(y2, ambient, threshold=0.5)
    expect_true(scaling > scaling2)

    scaling3 <- maximumAmbience(y2, ambient, threshold=0.01)
    expect_true(scaling < scaling3)

    scaling4 <- maximumAmbience(y2, ambient, dispersion=0.5)
    expect_true(scaling4 > scaling)
})

test_that("maximumAmbience function is robust to extra zeroes", {
    scaling <- maximumAmbience(y2, ambient)
    scaling2 <- maximumAmbience(c(0,0,0,y2), c(0,0,0,ambient))
    expect_equal(scaling, scaling2)

    prof <- maximumAmbience(y2, ambient, mode="profile")
    prof2 <- maximumAmbience(c(0,0,0,y2), c(0,0,0,ambient), mode="profile")
    expect_equal(c(0,0,0,prof), prof2[,1])

    prop <- maximumAmbience(y2, ambient, mode="proportion")
    prop2 <- maximumAmbience(c(0,0,0,y2), c(0,0,0,ambient), mode="proportion")
    expect_equal(c(NaN,NaN,NaN,prop), prop2[,1])
})

test_that("maximumAmbience function preserves names", {
    names(y2) <- seq_along(y2)
    expect_warning(maximumAmbience(y2, ambient, mode="profile"), "do not have the same feature names")

    names(ambient) <- seq_along(y2)
    prop <- maximumAmbience(y2, ambient, mode="profile")
    expect_identical(rownames(prop), names(y2))

    prop <- maximumAmbience(y2, ambient, mode="proportion")
    expect_identical(rownames(prop), names(y2))
    
    prop <- maximumAmbience(cbind(A=y2), ambient)
    expect_identical(names(prop), "A")
})
