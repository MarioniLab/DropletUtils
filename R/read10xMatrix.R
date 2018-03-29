#' @export
#' @importFrom Matrix readMM sparseMatrix
#' @importFrom methods as
#' @importFrom utils tail
# Note that this function is the reason for Suggests: HDF5Array
read10xMatrix <- function(file, hdf5.out=FALSE, chunk.size) {
    if (is.character(file)) { 
        fhandle <- file(file, open='r')
        on.exit(close(fhandle))
    } else {
        fhandle <- file
    }

    # Checking if we should chunk.
    if (!hdf5.out) { 
        out <- readMM(file)
        out <- as(out, "dgCMatrix")
        return(out)
    }

    chunk.size <- as.integer(chunk.size)
    if (chunk.size <= 0L) { 
        stop("chunk size should be a positive integer") 
    }

    # Checking the input type.
    type <- readLines(fhandle, 1)
    if (!grepl("%%MatrixMarket matrix coordinate (real|integer) general", type)) {
        stop("expected numeric/integer matrix in MatrixMarket format")
    }
    if (strsplit(type, split=" ")[[1]][4] == "integer") {
        input.x <- integer()
    } else {
        input.x <- double()
    }
    
    # Checking dimensions.
    dims <- scan(fhandle, what=integer(), nmax=3, comment.char = "%", quiet=TRUE)
    nr <- dims[1]
    nc <- dims[2]
    nz <- dims[3]

    out <- .Call(cxx_load_tenx_to_hdf5, fhandle, chunk.size, input.x, nr, nc, nz)
    return(out)
}

.scan_mm_file <- function(fhandle, chunk.size, x.type, order=FALSE) {
    out <- scan(fhandle, nmax = chunk.size, quiet = TRUE, what = list(i = integer(), j = integer(), x = x.type)) 
    if (order) {
        o <- order(out$j)
        out$i <- out$i[o]
        out$j <- out$j[o]
        out$x <- out$x[o]
    }
    return(out)
}

