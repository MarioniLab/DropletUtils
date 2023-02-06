#' @importFrom stats rexp rgamma runif rpois
#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
simCounts <- function(ngenes=100, nempty=10000, empty.prof=seq_len(ngenes), empty.rate=0.04,
                      nsmall=100, small.prof=runif(ngenes), small.shape=20, small.rate=0.1,
                      nlarge=1000, large.prof=empty.prof, large.shape=10, large.rate=0.01) 
# This simulates some counts for use in testing emptyDrops, and other 
# functions related to the count matrix (e.g., downsampleMatrix).
#     
# written by Aaron Lun
# created 6 January 2018
{ 
   
    # Normalizing the profiles. 
    empty.prof <- empty.prof/sum(empty.prof)
    large.prof <- large.prof/sum(large.prof)
    small.prof <- small.prof/sum(small.prof)

    # Simulating empty counts.
    total.count <- rexp(nempty, rate=empty.rate)
    empty.counts <- matrix(rpois(ngenes*nempty, lambda=outer(empty.prof, total.count)), ncol=nempty, nrow=ngenes)
    empty.counts <- as(empty.counts, "CsparseMatrix")

    # Simulating large counts.
    total.count <- rgamma(nlarge, shape=large.shape, rate=large.rate)
    large.counts <- matrix(rpois(ngenes*nlarge, lambda=outer(large.prof, total.count)), ncol=nlarge, nrow=ngenes)
    large.counts <- as(large.counts, "CsparseMatrix")

    # Simulating small counts.
    total.count <- rgamma(nsmall, shape=small.shape, rate=small.rate)
    small.counts <- matrix(rpois(ngenes*nsmall, lambda=outer(small.prof, total.count)), ncol=nsmall, nrow=ngenes)
    small.counts <- as(small.counts, "CsparseMatrix")

    # Creating the output matrix.
    out <- cbind(empty.counts, large.counts, small.counts)
    rownames(out) <- paste0("GENE", seq_len(ngenes))
    return(out)
}

