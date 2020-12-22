# This tests the control ambience function.
# library(testthat); library(DropletUtils); source("test-conamb.R")

# Making up some data.
set.seed(1100000)
ambient <- c(runif(900, 0, 0.1), runif(100))
y <- rpois(1000, ambient * 50)
y2 <- y + rpois(1000, 5) # adding some actual biology.

test_that("controlAmbience function works with a variety of inputs", {
    ref <- controlAmbience(y, ambient, 1:100)
    expect_equal(unname(ref), sum(y[1:100])/sum(ambient[1:100]))

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
})

test_that("controlAmbience function works with gene sets", {
    sets <- list(1:100, 200:300)
    ref <- controlAmbience(y, ambient, sets)
    expect_equal(unname(ref), 
        min(
            sum(y[1:100])/sum(ambient[1:100]),
            sum(y[200:300])/sum(ambient[200:300])
        )
    )

    # Same with matrices.
    ref <- controlAmbience(cbind(y, y2), cbind(ambient, ambient+1), sets)
    expect_equal(unname(ref), 
        c(
            unname(controlAmbience(y, ambient, sets)),
            unname(controlAmbience(y2, ambient+1, sets))
        )
    )
})

test_that("controlAmbience function works with a variety of outputs", {
    out <- controlAmbience(y, ambient, 1:100, mode="profile")
    expect_identical(ncol(out), 1L)
    expect_identical(nrow(out), length(y))

    out <- controlAmbience(y, ambient, 1:100, mode="proportion")
    expect_identical(ncol(out), 1L)
    expect_identical(nrow(out), length(y))

    # Behaves correctly with zeroes.
    has.zero <- controlAmbience(c(0, y), c(0, ambient), 1+1:100, mode="proportion")
    expect_identical(has.zero[1], NaN)
    expect_identical(has.zero[-1,,drop=FALSE], out)
})

test_that("controlAmbience respects names", {
    names(y2) <- seq_along(y2)
    expect_warning(controlAmbience(y2, ambient, 1:100, mode="profile"), "are not the same")
    expect_warning(controlAmbience(cbind(y2), cbind(ambient), 1:100, mode="profile"), "do not have the same row names")

    names(ambient) <- seq_along(y2)
    prop <- controlAmbience(y2, ambient, 1:100, mode="profile")
    expect_identical(rownames(prop), names(y2))

    prop <- controlAmbience(y2, ambient, 1:100, mode="proportion")
    expect_identical(rownames(prop), names(y2))
    
    prop <- controlAmbience(cbind(A=y2), ambient, 1:100)
    expect_identical(names(prop), "A")
})
