#' @export
write10xCounts <- function(path, x, barcodes=colnames(x), gene.id=rownames(x), gene.symbol=gene.id, gene.type="Gene Expression",
    overwrite=FALSE, type=c("auto", "sparse", "HDF5"), genome="unknown", version=c("2", "3"))
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
    version <- match.arg(version)
    type <- .type_chooser(path, match.arg(type))
    if (type=="sparse") {
        .write_sparse(temp.path, x, barcodes, gene.id, gene.symbol, gene.type, version=version)
    } else {
        .write_hdf5(temp.path, genome, x, barcodes, gene.id, gene.symbol, gene.type, version=version)
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
#' @importFrom R.utils gzip
.write_sparse <- function(path, x, barcodes, gene.id, gene.symbol, gene.type, version="2") {
    dir.create(path, showWarnings=FALSE)
    gene.info <- data.frame(gene.id, gene.symbol, stringsAsFactors=FALSE)

    if (version=="3") {
        gene.info$gene.type <- rep(gene.type, length.out=nrow(gene.info))
        mhandle <- file.path(path, "matrix.mtx")
        bhandle <- gzfile(file.path(path, "barcodes.tsv.gz"), open="wb")
        fhandle <- gzfile(file.path(path, "features.tsv.gz"), open="wb")
        on.exit({
            close(bhandle)
            close(fhandle)
        })
    } else {
        mhandle <- file.path(path, "matrix.mtx")
        bhandle <- file.path(path, "barcodes.tsv")
        fhandle <- file.path(path, "genes.tsv")
    }

    writeMM(x, file=mhandle)
    write(barcodes, file=bhandle)
    write.table(gene.info, file=fhandle, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

    if (version=="3") {
        # Annoyingly, writeMM doesn't take connection objects.
        gzip(mhandle)
    }

    return(NULL)
}

#' @importFrom rhdf5 h5createFile h5createGroup h5write
#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
.write_hdf5 <- function(path, genome, x, barcodes, gene.id, gene.symbol, gene.type, version="3") {
    h5createFile(path)

    if (version=="3") {
        group <- "matrix"
    } else {
        group <- genome
    }
    h5createGroup(path, group)

    h5write(barcodes, file=path, name=paste0(group, "/barcodes"))

    # Saving feature information.
    if (version=="3") {
        h5createGroup(path, file.path(group, "features"))
        h5write(gene.id, file=path, name=paste0(group, "/features/id"))
        h5write(gene.symbol, file=path, name=paste0(group, "/features/name"))

        h5write(rep(gene.type, length.out=length(gene.id)),
            file=path, name=paste0(group, "/features/feature_type"))

        h5write(rep(genome, length.out=length(gene.id)),
            file=path, name=paste0(group, "/features/genome"))

    } else {
        h5write(gene.id, file=path, name=paste0(group, "/genes"))
        h5write(gene.symbol, file=path, name=paste0(group, "/gene_names"))
    }

    # Saving matrix information.
    x <- as(x, "dgCMatrix")
    h5write(x@x, file=path, name=paste0(group, "/data"))
    h5write(dim(x), file=path, name=paste0(group, "/shape"))
    if (version=="3") {
        ind <- x@i - 1L # zero indexed in version 3.0.
    } else {
        ind <- x@i
    }
    h5write(ind, file=path, name=paste0(group, "/indices"))
    h5write(x@p, file=path, name=paste0(group, "/indptr"))

    return(NULL)
}

.type_chooser <- function(path, type) {
    if (type=="auto") {
        type <- if (grepl("\\.h5", path)) "HDF5" else "sparse"
    }
    type
}
