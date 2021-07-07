# This tests the control ambience function.
# library(testthat); library(DropletUtils); source("test-amb-control.R")

# Making up some data.
set.seed(1100000)
ambient <- c(runif(900, 0, 0.1), runif(100))
y <- rpois(1000, ambient * 50)
y2 <- y + rpois(1000, 5) # adding some actual biology.

test_that("ambientContribNegative function works with a variety of inputs", {
    ref <- ambientContribNegative(y, ambient, 1:100)
    expect_equal(unname(ref), sum(y[1:100])/sum(ambient[1:100]))

    out <- ambientContribNegative(cbind(y), ambient, 1:100)
    expect_identical(ref, out)
    expect_identical(length(ref), 1L)

    ref2 <- ambientContribNegative(y2, ambient, 1:100)
    out2 <- ambientContribNegative(cbind(y, y2), ambient, 1:100)
    expect_identical(ref, out2[1])
    expect_equivalent(ref2, out2[2])

    ref2x <- ambientContribNegative(y2, ambient+1, 1:100)
    out2x <- ambientContribNegative(cbind(y, y2), cbind(ambient, ambient+1), 1:100)
    expect_identical(ref, out2x[1])
    expect_equivalent(ref2x, out2x[2])
})

test_that("ambientContribNegative function works with gene sets", {
    sets <- list(1:100, 200:300)
    ref <- ambientContribNegative(y, ambient, sets)
    expect_equal(unname(ref), 
        min(
            sum(y[1:100])/sum(ambient[1:100]),
            sum(y[200:300])/sum(ambient[200:300])
        )
    )

    # Same with matrices.
    ref <- ambientContribNegative(cbind(y, y2), cbind(ambient, ambient+1), sets)
    expect_equal(unname(ref), 
        c(
            unname(ambientContribNegative(y, ambient, sets)),
            unname(ambientContribNegative(y2, ambient+1, sets))
        )
    )

    # Works for more than one gene set.
    sets <- list(1:100, 200:300, 400:500)
    ref <- ambientContribNegative(y, ambient, sets)
    expect_equal(unname(ref), 
        median(
            c(sum(y[1:100])/sum(ambient[1:100]),
            sum(y[200:300])/sum(ambient[200:300]),
            sum(y[400:500])/sum(ambient[400:500]))
        )
    )

})

test_that("ambientContribNegative function works with a variety of outputs", {
    out <- ambientContribNegative(y, ambient, 1:100, mode="profile")
    expect_identical(ncol(out), 1L)
    expect_identical(nrow(out), length(y))

    out <- ambientContribNegative(y, ambient, 1:100, mode="proportion")
    expect_identical(ncol(out), 1L)
    expect_identical(nrow(out), length(y))

    # Behaves correctly with zeroes.
    has.zero <- ambientContribNegative(c(0, y), c(0, ambient), 1+1:100, mode="proportion")
    expect_identical(has.zero[1], NaN)
    expect_identical(has.zero[-1,,drop=FALSE], out)
})

test_that("ambientContribNegative respects names", {
    names(y2) <- seq_along(y2)
    expect_warning(ambientContribNegative(y2, ambient, 1:100, mode="profile"), "are not the same")
    expect_warning(ambientContribNegative(cbind(y2), cbind(ambient), 1:100, mode="profile"), "do not have the same row names")

    names(ambient) <- seq_along(y2)
    prop <- ambientContribNegative(y2, ambient, 1:100, mode="profile")
    expect_identical(rownames(prop), names(y2))

    prop <- ambientContribNegative(y2, ambient, 1:100, mode="proportion")
    expect_identical(rownames(prop), names(y2))
    
    prop <- ambientContribNegative(cbind(A=y2), ambient, 1:100)
    expect_identical(names(prop), "A")
})
