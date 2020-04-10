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
    expect_equal(c(0,0,0,prof), prof2)

    prop <- maximumAmbience(y2, ambient, mode="proportion")
    prop2 <- maximumAmbience(c(0,0,0,y2), c(0,0,0,ambient), mode="proportion")
    expect_equal(c(NaN,NaN,NaN,prop), prop2)
})
