#' @export
#' @importFrom utils write.table
#' @importFrom Matrix writeMM
write10xCounts <- function(path, x, barcodes=colnames(x), gene.id=rownames(x), gene.symbol=gene.id, overwrite=FALSE)
# Writes a count matrix to the specified path in the 10x style.   
# This allows us to easily create things for testing read10xCounts. 
# 
# written by Aaron Lun
# created 6 January 2018
{
    # Creating a temp folder at the local path, as we have permissions there.
    # This avoids problems with 'path' already existing.
    temp.path <- tempfile(tmpdir=dirname(path)) 
    dir.create(temp.path, showWarnings=FALSE)
    on.exit({ 
        if (file.exists(temp.path)) { unlink(temp.path, recursive=TRUE) } 
    })

    # Saving all of the various bits and pieces.
    writeMM(x, file=file.path(temp.path, "matrix.mtx"))

    if (ncol(x)!=length(barcodes)) { 
        stop("'barcodes' must of of the same length as 'ncol(x)'")
    }
    write(barcodes, file=file.path(temp.path, "barcodes.tsv"))

    if (length(gene.id)!=length(gene.symbol) || length(gene.id)!=nrow(x)) {
        stop("lengths of 'gene.id' and 'gene.symbol' must be equal to 'nrow(x)'")
    }
    write.table(data.frame(gene.id, gene.symbol), file=file.path(temp.path, "genes.tsv"),
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")    

    # We don't put this at the top as the write functions might fail; 
    # in which case, we would have deleted the existing 'path' for nothing.
    if (overwrite) {
        unlink(path, recursive=TRUE)
    } else if (file.exists(path)) { 
        stop("specified 'path' already exists")
    }
    file.rename(temp.path, path)
    return(invisible(TRUE))
}
