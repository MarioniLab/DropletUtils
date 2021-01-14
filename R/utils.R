.rounded_to_integer <- function(m, round=TRUE) {
    if (round) {
        m <- round(m)
    }
    m
}

#' @importFrom Matrix colSums
.intColSums <- function(m, BPPARAM) {
    # Enforcing discreteness mainly for emptyDrops()'s Monte Carlo step.
    as.integer(round(colSums(m)))
}

#' @importFrom DelayedArray DelayedArray SparseArraySeed which
#' @importClassesFrom DelayedArray DelayedArray SparseArraySeed
.realize_da_to_memory <- function(m) {
    if (is(m, "DelayedArray")) {
        if (!is(seed(m), "SparseArraySeed")) {
            idx <- which(m!=0, arr.ind=TRUE)
            m <- DelayedArray(
                SparseArraySeed(
                    nzindex=idx, 
                    nzdata=m[idx],
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

