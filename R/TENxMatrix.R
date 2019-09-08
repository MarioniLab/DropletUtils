################################
# Defining the TENxMatrixSeed class.

#' @export
#' @importClassesFrom DelayedArray Array
setClass("TENxMatrixSeed",
    contains="Array",
    representation(
        filepath="character",    # Absolute path to the HDF5 file so the
                                 # object doesn't break when the user
                                 # changes the working directory (e.g. with
                                 # setwd()).
        group="character",       # Name of the group in the HDF5 file
                                 # containing the 10x Genomics data.
        dim="integer",
        dimnames="list",
        col_ranges="data.frame"  # Can't use an IRanges object for this at the
                                 # moment because they don't support Linteger
                                 # start/end values yet.
    )
)

################################
# Getters that do not touch the file.

#' @export
#' @importFrom DelayedArray path
setMethod("path", "TENxMatrixSeed", function(object) object@filepath)

#' @export
setMethod("dim", "TENxMatrixSeed", function(x) x@dim)

#' @export
setMethod("dimnames", "TENxMatrixSeed", function(x) {
    if (is.null(x@dimnames[[1]]) && is.null(x@dimnames[[2]])) {
        NULL
    } else {
        x@dimnames
    }
})

#' @export
#' @importFrom DelayedArray chunkdim
setMethod("chunkdim", "TENxMatrixSeed", function(x) 
### Chunks are defined column-wide, consistent with a column-sparse
### compressed matrix.
{ 
    c(nrow(x), 1L)
})

################################
# Getters that pull from disk.

#' @importFrom rhdf5 h5read
.get_TENx_component <- function(filepath, group, name, idx=NULL) {
    name <- paste0(group, "/", name)
    if (!is.null(idx))
        idx <- list(idx)
    as.vector(h5read(filepath, name, index=idx))
}

.get_shape <- function(filepath, group) {
    .get_TENx_component(filepath, group, "shape")
}

#' @importFrom rhdf5 h5read
.get_indptr <- function(filepath, group) {
    name <- paste0(group, "/indptr")
    as.vector(h5read(filepath, name, bit64conversion="double"))
}

.get_data <- function(filepath, group, idx=NULL) {
    .get_TENx_component(filepath, group, "data", idx=idx)
}

#' @importFrom rhdf5 h5read
.linear_get_data <- function(filepath, group, start=NULL, count=NULL) {
    name <- paste0(group, "/data")
    as.vector(h5read(filepath, name, start=start, count=count))
}

.get_row_indices <- function(filepath, group, idx=NULL) .get_TENx_component(filepath, group, "indices", idx=idx)

#' @importFrom rhdf5 h5read
.linear_get_row_indices <- function(filepath, group, start=NULL, count=NULL) {
    name <- paste0(group, "/indices")
    as.vector(h5read(filepath, name, start=start, count=count))
}

.get_genes <- function(filepath, group, idx=NULL) {
    if (!HDF5Array:::h5exists(filepath, paste0(group, "/genes"))) {
        NULL
    } else {
        .get_TENx_component(filepath, group, "genes", idx=idx)
    }
}

.get_barcodes <- function(filepath, group, idx=NULL) {
    if (!HDF5Array:::h5exists(filepath, paste0(group, "/barcodes"))) {
        NULL
    } else {
        .get_TENx_component(filepath, group, "barcodes", idx=idx)
    }
}

################################

#' @importFrom S4Vectors normalizeSingleBracketSubscript
#' @importFrom utils relist
#' @importFrom IRanges PartitioningByWidth
.get_data_indices_by_col <- function(x, j)
# 'j' must be an integer vector containing valid col indices.
# Return data indices in a NumericList object parallel to 'j' i.e. with
# one list element per col index in 'j'.
{
    j <- normalizeSingleBracketSubscript(j, x, exact = FALSE)
    offset <- col_ranges[j, "start"] - 1L
    lengths <- col_ranges[j, "width"]

    lengths_len <- length(lengths)
    if (lengths_len == 0L) {
        idx2 <- numeric(0)
    } else {
        offsets <- offset - cumsum(c(0L, lengths[-lengths_len]))
        idx2 <- seq_len(sum(lengths)) + rep.int(offsets, lengths)
    }

    relist(idx2, PartitioningByWidth(lengths))
}

################################

#' @importFrom DelayedArray SparseArraySeed
#' @importFrom utils relist
#' @importFrom IRanges PartitioningByWidth
.extract_data_from_adjacent_cols <- function(x, j1, j2, as.sparse=FALSE)
### 'j1' and 'j2' must be 2 integers representing a valid range of col indices.
### If 'as.sparse=FALSE', return a NumericList or IntegerList object parallel
### to 'j1:j2' i.e. with one list element per col index in 'j1:j2'.
### If 'as.sparse=TRUE', return a SparseArraySeed object.
{
    start <- x@col_ranges[j1, "start"]

    j12 <- j1:j2
    count_per_col <- x@col_ranges[j12, "width"]
    count <- sum(count_per_col)
    ans_nzdata <- .linear_get_data(x@filepath, x@group, start=start, count=count)

    if (!as.sparse) {
        return(relist(ans_nzdata, PartitioningByWidth(count_per_col)))
    }

    row_indices <- .linear_get_row_indices(x@filepath, x@group, start=start, count=count) + 1L
    col_indices <- rep.int(j12, count_per_col)
    ans_aind <- cbind(row_indices, col_indices, deparse.level=0L)
    SparseArraySeed(dim(x), ans_aind, ans_nzdata, check=FALSE)
}

################################

#' @importFrom utils relist
.random_extract_nonzero_data_by_col <- function(x, j)
### Extract nonzero data using the "random" method. This method is
### based on h5read( , index=idx) which retrieves an arbitrary/random
### subset of the data.
### 'j' must be an integer vector containing valid col indices. It cannot
### be NULL.
{
    data_indices <- .get_data_indices_by_col(x, j)
    idx2 <- unlist(data_indices, use.names=FALSE)
    data <- .get_data(x@filepath, x@group, idx=idx2)
    relist(data, data_indices)
}

.linear_extract_nonzero_data_by_col <- function(x, j)
### Extract nonzero data using the "linear" method. This method is
### based on h5read( , start=start, count=count) which retrieves a
### linear subset of the data and is much faster than doing
### h5read( , index=list(seq(start, length.out=count))).
### 'j' must be NULL or an integer vector containing valid col indices. It
### should not be empty.
{
    if (is.null(j)) {
        j1 <- 1L
        j2 <- ncol(x)
    } else {
        stopifnot(is.numeric(j), length(j) != 0L)
        j1 <- min(j)
        j2 <- max(j)
    }
    nonzero_data <- .extract_data_from_adjacent_cols(x, j1, j2)
    if (is.null(j)) { 
        nonzero_data
    } else {
        nonzero_data[j - j1 + 1L]
    }
}

.normarg_method <- function(method, j) {
    if (method != "auto") {
        return(method)
    }

    if (is.null(j)) {
        return("linear")
    }

    if (length(j) == 0L) {
        return("random")
    }

    j1 <- min(j)
    j2 <- max(j)

    ## 'ratio' is > 0 and <= 1. A value close to 1 indicates that the columns
    ## to extract are close from each other (a value of 1 indicating that
    ## they are adjacent e.g. j <- 18:25). A value close to 0 indicates that
    ## they are far apart from each other i.e. that they are separated by many
    ## columns that are not requested. The "linear" method is very efficient
    ## when 'ratio' is close to 1. It is so much more efficient than the
    ## "random" method (typically 10x or 20x faster) that we choose it when
    ## 'ratio' is >= 0.2
    ratio <- length(j) / (j2 - j1 + 1L)
    if (ratio >= 0.2) { 
        "linear" 
    } else {
        "random"
    }
}

.extract_nonzero_data_by_col <- function(x, j, method=c("auto", "random", "linear"))
### 'j' must be NULL or an integer vector containing valid col indices.
### Return a NumericList or IntegerList object parallel to 'j' i.e. with
### one list element per col index in 'j'.
{
    method <- match.arg(method)
    method <- .normarg_method(method, j)
    if (method == "random") {
        .random_extract_nonzero_data_by_col(x, j)
    } else {
        .linear_extract_nonzero_data_by_col(x, j)
    }
}

################################

#' @importFrom DelayedArray SparseArraySeed
.random_load_SparseArraySeed_from_TENxMatrixSeed <- function(x, i, j)
### Load sparse data using the "random" method.
### This method is based on h5read( , index=idx) which retrieves an
### arbitrary/random subset of the data.
### 'i' must be NULL or an integer vector containing valid row indices.
### 'j' must be an integer vector containing valid col indices. It cannot
### be NULL.
### Both 'i' and 'j' can contain duplicates. Duplicates in 'i' have no effect
### on the output but duplicates in 'j' will produce duplicates in the output.
### Return a SparseArraySeed object.
{
    stopifnot(is.null(i) || is.numeric(i), is.numeric(j))

    data_indices <- .get_data_indices_by_col(x, j)
    idx2 <- unlist(data_indices, use.names=FALSE)
    row_indices <- .get_row_indices(x@filepath, x@group, idx=idx2) + 1L
    col_indices <- rep.int(j, lengths(data_indices))

    if (!is.null(i)) {
        keep_idx <- which(row_indices %in% i)
        idx2 <- idx2[keep_idx]
        row_indices <- row_indices[keep_idx]
        col_indices <- col_indices[keep_idx]
    }

    ans_aind <- cbind(row_indices, col_indices, deparse.level=0L)
    ans_nzdata <- .get_data(x@filepath, x@group, idx=idx2)
    SparseArraySeed(dim(x), ans_aind, ans_nzdata, check=FALSE)
}

.linear_load_SparseArraySeed_from_TENxMatrixSeed <- function(x, j)
### Load sparse data using the "linear" method.
### This method is based on h5read( , start=start, count=count) which
### retrieves a linear subset of the data and is much faster than doing
### h5read( , index=list(seq(start, length.out=count))).
### 'j' must be NULL or a non-empty integer vector containing valid
### col indices.  ### The output is not affected by duplicates in 'j'.
### Return a SparseArraySeed object.
{
    if (is.null(j)) {
        j1 <- 1L
        j2 <- ncol(x)
    } else {
        stopifnot(is.numeric(j), length(j) != 0L)
        j1 <- min(j)
        j2 <- max(j)
    }
    .extract_data_from_adjacent_cols(x, j1, j2, as.sparse=TRUE)
}

.load_SparseArraySeed_from_TENxMatrixSeed <- function(x, index, method=c("auto", "random", "linear"))
### Duplicates in 'index[[1]]' are ok and won't affect the output.
### Duplicates in 'index[[2]]' are ok but might introduce duplicates
### in the output so should be avoided.
### Return a SparseArraySeed object.
{
    i <- index[[1L]]
    j <- index[[2L]]
    method <- match.arg(method)
    method <- .normarg_method(method, j)
    if (method == "random") {
        .random_load_SparseArraySeed_from_TENxMatrixSeed(x, i, j)
    } else {
        .linear_load_SparseArraySeed_from_TENxMatrixSeed(x, j)
    }
}

################################
# TENxMatrixSeed methods.

#' @export
#' @importFrom DelayedArray extract_array
setMethod("extract_array", "TENxMatrixSeed", function(x, index) {
    if (length(index[[1]]) == 0L || length(index[[2]]) == 0L) {
        data0 <- .get_data(x@filepath, x@group, idx=integer(0))
        array(data0, dim=ans_dim)
    } else {
        sas <- .load_SparseArraySeed_from_TENxMatrixSeed(x, index) 
        extract_array(sas, index)
    }
})

#' @export
#' @importFrom S4Vectors wmsg isSingleString
TENxMatrixSeed <- function(filepath, group)
{
    filepath <- HDF5Array:::normarg_path(filepath, "'filepath'", "10x Genomics dataset")
    if (!isSingleString(group)) {
        stop(wmsg("'group' must be a single string specifying the name ",
                  "of the group in the HDF5 file containing the ",
                  "10x Genomics data"))
    }
    if (group == "") {
        stop(wmsg("'group' cannot be the empty string"))
    }

    ## dim
    dim <- .get_shape(filepath, group)
    stopifnot(length(dim) == 2L)

    ## dimnames
    rownames <- .get_genes(filepath, group)
    stopifnot(is.null(rownames) || length(rownames) == dim[[1L]])
    colnames <- .get_barcodes(filepath, group)
    stopifnot(is.null(colnames) || length(colnames) == dim[[2L]])
    dimnames <- list(rownames, colnames)

    ## col_ranges
    data_len <- HDF5Array:::h5length(filepath, paste0(group, "/data"))
    indices_len <- HDF5Array:::h5length(filepath, paste0(group, "/indices"))
    stopifnot(data_len == indices_len)

    indptr <- .get_indptr(filepath, group)
    stopifnot(length(indptr) == dim[[2L]] + 1L,
        indptr[[0L]] == 0L,
        indptr[[length(indptr)]] == data_len)

    col_ranges <- data.frame(start=indptr[-length(indptr)] + 1,
        width=as.integer(diff(indptr)))

    new2("TENxMatrixSeed", filepath=filepath,
        group=group,
        dim=dim,
        dimnames=dimnames,
        col_ranges=col_ranges)
}


setMethod("show", "TENxMatrixSeed", function(object) {
    cat(DelayedArray:::array_as_one_line_summary(object), ":\n", sep="")
    cat("# dirname: ", dirname(object), "\n", sep="")
    cat("# basename: ", basename(object), "\n", sep="")
    cat("# group: ", object@group, "\n", sep="")
})

################################
# Sparsity-related methods.

#' @export
#' @importFrom DelayedArray sparsity
setMethod("sparsity", "TENxMatrixSeed", function(x) {
    data_len <- HDF5Array:::h5length(x@filepath, paste0(x@group, "/data"))
    1 - data_len / length(x)
})

#' @export
#' @importFrom DelayedArray is_sparse
setMethod("is_sparse", "TENxMatrixSeed", function(x) TRUE)

#' @export
#' @importFrom DelayedArray extract_sparse_array
setMethod("extract_sparse_array", "TENxMatrixSeed", function(x, index)
{
    sas <- .load_SparseArraySeed_from_TENxMatrixSeed(x, index)
    extract_sparse_array(sas, index)
})

#' @export
#' @importFrom DelayedArray read_sparse_block
setMethod("read_sparse_block", "TENxMatrixSeed", function(x, viewport)
### The default "read_sparse_block" method defined in DelayedArray would work
### just fine on a TENxMatrixSeed object (thanks to the "extract_sparse_array"
### method for TENxMatrixSeed objects defined above), but we overwrite it with
### the method below which is slightly more efficient. That's because the
### method below calls read_sparse_block() on the SparseArraySeed object
### returned by .load_SparseArraySeed_from_TENxMatrixSeed() and this is
### faster than calling extract_sparse_array() on the same object (which
### is what the "extract_sparse_array" method for TENxMatrixSeed would do
### when called by the default "read_sparse_block" method).
### Not sure the difference is significant enough for this extra method to
### be worth it though, because, time is really dominated by I/O here, that
### is, by the call to .load_SparseArraySeed_from_TENxMatrixSeed().
{
    index <- makeNindexFromArrayViewport(viewport, expand.RangeNSBS=TRUE)
    sas <- .load_SparseArraySeed_from_TENxMatrixSeed(x, index)  # I/O

    ## Unlike the "extract_sparse_array" method for TENxMatrixSeed objects
    ## defined above, we use read_sparse_block() here, which is faster than
    ## using extract_sparse_array().
    read_sparse_block(sas, viewport)  # in-memory
})

################################
# Coercion to dgCMatrix.

#' @importFrom Matrix sparseMatrix
.from_TENxMatrixSeed_to_dgCMatrix <- function(from) {
    row_indices <- .get_row_indices(from@filepath, from@group) + 1L
    indptr <- .get_indptr(from@filepath, from@group)
    data <- .get_data(from@filepath, from@group)
    sparseMatrix(i=row_indices, p=indptr, x=data, dims=dim(from),
        dimnames=dimnames(from))
}

#' @importClassesFrom Matrix dgCMatrix
setAs("TENxMatrixSeed", "dgCMatrix", .from_TENxMatrixSeed_to_dgCMatrix)

#' @importClassesFrom Matrix sparseMatrix
setAs("TENxMatrixSeed", "sparseMatrix", .from_TENxMatrixSeed_to_dgCMatrix)

################################
# Defining the TENxMatrix class. This is done mostly for cosmetic reasons, to
# hide the DelayedMatrix class from the user. Thus, the user will see and 
# manipulate TENxMatrix objects instead of DelayedMatrix objects.

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setClass("TENxMatrix",
    contains="DelayedMatrix",
    representation(seed="TENxMatrixSeed")
)

#' @export
#' @importFrom DelayedArray DelayedArray new_DelayedArray
setMethod("DelayedArray", "TENxMatrixSeed",
    function(seed) new_DelayedArray(seed, Class="TENxMatrix")
)

#' @export
#' @importFrom DelayedArray DelayedArray
TENxMatrix <- function(filepath, group) {
    if (is(filepath, "TENxMatrixSeed")) {
        if (!missing(group)) {
            stop(wmsg("TENxMatrix() must be called with a single argument ",
                      "when passed a TENxMatrixSeed object"))
        }
        seed <- filepath
    } else {
        seed <- TENxMatrixSeed(filepath, group)
    }
    DelayedArray(seed)
}

#' @export
#' @importFrom DelayedArray sparsity
setMethod("sparsity", "TENxMatrix", function(x) sparsity(x@seed))

#' @export
#' @importFrom DelayedArray read_sparse_block
setMethod("read_sparse_block", "TENxMatrix", function(x, viewport) {
    read_sparse_block(x@seed, viewport)
})

#' @importClassesFrom Matrix dgCMatrix
.from_TENxMatrix_to_dgCMatrix <- function(from) as(from@seed, "dgCMatrix")

#' @importClassesFrom Matrix dgCMatrix
setAs("TENxMatrix", "dgCMatrix", .from_TENxMatrix_to_dgCMatrix)

#' @importClassesFrom Matrix sparseMatrix
setAs("TENxMatrix", "sparseMatrix", .from_TENxMatrix_to_dgCMatrix)
