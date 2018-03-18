#' @export
#' @importFrom Matrix readMM sparseMatrix
#' @importFrom methods as
#' @importFrom utils tail
# Note that this function is the reason for Suggests: HDF5Array
read10xMatrix <- function(file, chunk.size=NA, hdf5.out=FALSE) {
    if (is.character(file)) { 
        fhandle <- file(file, open='r')
        on.exit(close(fhandle))
    } else {
        fhandle <- file
    }

    # Checking the input type.
    type <- readLines(fhandle, 1)
    if (type!="%%MatrixMarket matrix coordinate integer general") {
        stop("expected integer matrix in MatrixMarket format")
    }

    # Checking if we should chunk.
    if (is.na(chunk.size)) { 
        out <- readMM(file)
        out <- as(out, "dgCMatrix")
        return(out)
    }

    chunk.size <- as.integer(chunk.size)
    if (chunk.size <= 0L) { 
        stop("chunk size should be a positive integer") 
    }
    
    # Checking dimensions.
    dims <- scan(fhandle, what=integer(), nmax=3, comment.char = "%", quiet=TRUE)
    nr <- dims[1]
    nc <- dims[2]
    nz <- dims[3]

    if (!hdf5.out) { 
        collected.i <- collected.j <- integer(nz)
        collected.x <- numeric(nz)
        counter <- 0L

        while (counter < nz) {
            out <- .scan_mm_file(fhandle, chunk.size=chunk.size)
            idx <- (counter+1L):min(nz, counter+chunk.size)
            collected.i[idx] <- out$i
            collected.j[idx] <- out$j
            collected.x[idx] <- out$x
            counter <- tail(idx, 1L)
        }

        out <- sparseMatrix(i=collected.i, j=collected.j, x=collected.x, dims=c(nr, nc), giveCsparse=TRUE)
    } else {
        out <- .Call(cxx_load_tenx_to_hdf5, fhandle, chunk.size, nr, nc, nz)
    }
        
    return(out)
}

.scan_mm_file <- function(fhandle, chunk.size, order=FALSE) {
    out <- scan(fhandle, nmax = chunk.size, quiet = TRUE, what = list(i = integer(), j = integer(), x = integer()))
    if (order) {
        o <- order(out$j, out$i)
        out$i <- out$i[o]
        out$j <- out$j[o]
        out$x <- out$x[o]
    }
    return(out)
}

