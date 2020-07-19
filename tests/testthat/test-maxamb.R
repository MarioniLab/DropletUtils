# This tests the maximum ambience function.
# library(testthat); library(DropletUtils); source("test-maxamb.R")

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
    expect_identical(stuff[1], maximumAmbience(y, ambient))
    expect_identical(stuff[2], maximumAmbience(y2, ambient))

    ambient2 <- c(runif(900, 0, 0.1), runif(100))
    stuff <- maximumAmbience(cbind(y, y2), cbind(ambient, ambient2))
    expect_identical(stuff[1], maximumAmbience(y, ambient))
    expect_identical(stuff[2], maximumAmbience(y2, ambient2))

    stuff <- maximumAmbience(cbind(y, y2), cbind(ambient, ambient2), mode="profile")
    expect_identical(stuff[,1,drop=FALSE], maximumAmbience(y, ambient, mode="profile"))
    expect_identical(stuff[,2,drop=FALSE], maximumAmbience(y2, ambient2, mode="profile"))
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

test_that("controlAmbience function works with a variety of inputs and outputs", {
    ref <- controlAmbience(y, ambient, 1:100)
    out <- controlAmbience(cbind(y), ambient, 1:100)
    expect_identical(ref, out)
    expect_identical(length(ref), 1L)

    ref2 <- controlAmbience(y2, ambient, 1:100)
    out2 <- controlAmbience(cbind(y, y2), ambient, 1:100)
    expect_identical(ref, out2[1])
    expect_equivalent(ref2, out2[2])

    ref2x <- controlAmbience(y2, ambient+1, 1:100)
    out2x <- controlAmbience(cbind(y, y2), cbind(ambient, ambient+1), 1:100)
    expect_identical(ref, out2x[1])
    expect_equivalent(ref2x, out2x[2])

    # Checking the other output formats.
    out <- controlAmbience(y, ambient, 1:100, mode="profile")
    expect_identical(ncol(out), 1L)
    expect_identical(nrow(out), length(y))

    out <- controlAmbience(y, ambient, 1:100, mode="proportion")
    expect_identical(ncol(out), 1L)
    expect_identical(nrow(out), length(y))

    has.zero <- controlAmbience(c(0, y), c(0, ambient), 1+1:100, mode="proportion")
    expect_identical(has.zero[1], NaN)
    expect_identical(has.zero[-1], out)
})

