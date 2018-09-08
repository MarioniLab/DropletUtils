#' @export
write10xCounts <- function(path, x, barcodes=colnames(x), gene.id=rownames(x), gene.symbol=gene.id, 
    overwrite=FALSE, type=c("auto", "sparse", "HDF5"), group="group")
# Writes a count matrix to the specified path in the 10x style.   
# This allows us to easily create things for testing read10xCounts. 
# 
# written by Aaron Lun
# created 6 January 2018
{
    # Doing all the work on a temporary location next to 'path', as we have permissions there.
    # This avoids problems with 'path' already existing.
    temp.path <- tempfile(tmpdir=dirname(path)) 
    on.exit({ 
        if (file.exists(temp.path)) { unlink(temp.path, recursive=TRUE) } 
    })

    # Checking the values.
    if (length(gene.id)!=length(gene.symbol) || length(gene.id)!=nrow(x)) {
        stop("lengths of 'gene.id' and 'gene.symbol' must be equal to 'nrow(x)'")
    }
    if (ncol(x)!=length(barcodes)) { 
        stop("'barcodes' must of of the same length as 'ncol(x)'")
    }

    # Determining what format to save in.
    type <- .type_chooser(path, match.arg(type))
    if (type=="sparse") {
        .write_sparse(temp.path, x, barcodes, gene.id, gene.symbol)
    } else {
        .write_hdf5(temp.path, group, x, barcodes, gene.id, gene.symbol)
    }

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

#' @importFrom utils write.table
#' @importFrom Matrix writeMM
.write_sparse <- function(path, x, barcodes, gene.id, gene.symbol) {
    dir.create(path, showWarnings=FALSE)

    # Saving all of the various bits and pieces.
    writeMM(x, file=file.path(path, "matrix.mtx"))

    write(barcodes, file=file.path(path, "barcodes.tsv"))

    write.table(data.frame(gene.id, gene.symbol), file=file.path(path, "genes.tsv"),
        row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

    return(NULL)
}

#' @importFrom rhdf5 h5createFile h5createGroup h5write
#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
.write_hdf5 <- function(path, group, x, barcodes, gene.id, gene.symbol) {
    h5createFile(path)
    h5createGroup(path, group)
    h5write(barcodes, file=path, name=paste0(group, "/barcodes"))
    h5write(gene.id, file=path, name=paste0(group, "/genes"))
    h5write(gene.symbol, file=path, name=paste0(group, "/gene_names"))
   
    x <- as(x, "dgCMatrix")
    h5write(x@x, file=path, name=paste0(group, "/data"))
    h5write(x@i, file=path, name=paste0(group, "/indices"))
    h5write(x@p, file=path, name=paste0(group, "/indptr"))
    h5write(dim(x), file=path, name=paste0(group, "/shape"))

    return(NULL)
}

.type_chooser <- function(path, type) {
    if (type=="auto") {
        type <- if (grepl("\\.h5", path)) "HDF5" else "sparse"
    }
    type
}
