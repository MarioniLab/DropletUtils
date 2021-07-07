#' @export
#' @importFrom Matrix sparseMatrix
makeCountMatrix <- function(gene, cell, all.genes=NULL, all.cells=NULL, value=NULL) 
# This makes a count matrix from per-molecule information,
# namely the gene and cell to which each molecule is assigned.
#
# written by Aaron Lun
# created 21 December 2017    
{
    # Checking input lengths.
    stopifnot(identical(length(cell), length(gene)))
    if (is.null(value)) {
        value <- rep(1, length(gene))
    } else {
        stopifnot(identical(length(cell), length(value)))
    }
    
    # Checking the input gene.
    gene.out <- .check_indices(gene, all.genes, "gene")
    gene <- gene.out$val
    all.genes <- gene.out$tab
    ngenes <- gene.out$N

    # Checking the input cell.
    cell.out <- .check_indices(cell, all.cells, "cell")
    cell <- cell.out$val
    all.cells <- cell.out$tab
    ncells <- cell.out$N

    # Creating the matrix.
    out <- sparseMatrix(i=gene, j=cell, x=value, 
                        dims=c(ngenes, ncells),
                        use.last.ij=FALSE, # Adds up the duplicates.
                        giveCsparse=TRUE)
    rownames(out) <- all.genes
    colnames(out) <- all.cells
    return(out)
}

.check_indices <- function(val, tab, err) {
    if (is.character(val)) { 
        if (is.null(tab)) { 
            tab <- sort(unique(val))
        }
        val <- match(val, tab)
        if (any(is.na(val))) {
            stop(sprintf("entry of '%s' not in 'all.%ss'", err, err))
        }
        nvals <- length(tab)
    } else {
        if (!is.null(tab)) {
            nvals <- length(tab)
            if (length(val) && length(tab)<max(val)) {
                stop(sprintf("length of 'all.%ss' is less than 'max(%s)'", err, err))
            }
        } else {
            if (length(val)) { 
                nvals <- max(val)
            } else {
                nvals <- 0L
            }
        }
        if (length(val) && min(val)<=0) {
            stop(sprintf("indices in '%s' must be positive", err))
        }
    }
    return(list(val=val, tab=tab, N=nvals))
}
