.rounded_to_integer <- function(m, round=TRUE) {
    if (round) {
        m <- round(m)
    }
    m
}

#' @importFrom Matrix colSums
.intColSums <- function(m) {
    # Enforcing discreteness mainly for emptyDrops()'s Monte Carlo step.
    as.integer(round(colSums(m)))
}

#' @importFrom DelayedArray DelayedArray SparseArraySeed which seed
#' @importFrom beachmat whichNonZero
#' @importClassesFrom DelayedArray DelayedArray SparseArraySeed
.realize_DA_to_memory <- function(m, BPPARAM) {
    if (is(m, "DelayedArray")) {
        if (!is(seed(m), "SparseArraySeed")) {
            idx <- whichNonZero(m, BPPARAM)
            m <- DelayedArray(
                SparseArraySeed(
                    nzindex=cbind(idx$i, idx$j), 
                    nzdata=idx$x,
                    dim=dim(m),
                    dimnames=dimnames(m)
                )
            )
        }
    }

    m
}

#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
.parallelize <- function(BPPARAM) {
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    old
}
