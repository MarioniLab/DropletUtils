#' @export
#' @importFrom stats quantile
#' @importFrom utils head
defaultDrops <- function(m, expected = 3000, upper.quant = 0.99, lower.prop = 0.1)
# A function to call cells on library size, as performed by CellRanger
# 
# written by Jonathan Griffiths
# created 4 January 2018
{
    if(upper.quant > 1 | upper.quant < 0){
        stop("'upper.quant' should be a numeric value between 0 and 1")
    }
    
    if(lower.prop > 1 | lower.prop < 0){
        stop("'lower.prop' should be a numeric value between 0 and 1")
    }
    
    libs <- colSums(m)
    o <- order(libs, decreasing = TRUE)
    top <- libs[head(o, n = expected)]
    
    threshold <- quantile(top, upper.quant)*lower.prop
    return(libs > threshold)
}
