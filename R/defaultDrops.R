defaultDrops <- function(m, exp.cells = 3000, upper.pct = 0.99, lower.prop = 0.1)
# A function to call cells on library size, as performed by CellRanger
# 
# written by Jonathan Griffiths
# created 4 January 2018
{
  
  if(upper > 1 | upper < 0){
    stop("'upper' should be a numeric value between 0 and 1")
  }
  
  if(lower > 1 | lower < 0){
    stop("'lower' should be a numeric value between 0 and 1")
  }
  
  libs <- colSums(m)
  
  top <- libs[order(libs, decreasing = TRUE)[seq_len(exp.cells)]]
  
  threshold <- quantile(top, upper)*lower
  
  return(libs > threshold)
}
