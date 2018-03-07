#' @export
downsampleMatrix <- function(x, prop, bycol=TRUE) 
# Downsamples the count matrix to the specified proportion.
# This is done on a per-cell basis, possibly with cell-specific 'prop',
# unless 'bycol=FALSE', in which case it is done for the entire matrix.
#
# written by Aaron Lun
# created 18 December 2017    
{
    if (bycol) {
        prop <- rep(prop, length.out = ncol(x))
    }
    out <- .Call(cxx_downsample_matrix, x, prop, bycol)
    dimnames(out) <- dimnames(x)
    return(out)
}
