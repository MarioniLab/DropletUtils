downsampleCounts <- function(x, prop) 
# Downsamples the count matrix to the specified proportion
# (which can also be count specific).
#
# written by Aaron Lun
# created 18 December 2017    
{
    prop <- rep(prop, length.out = ncol(x))
    out <- .Call(cxx_downsample_matrix, x, prop)
    dimnames(out) <- dimnames(x)
    return(out)
}
