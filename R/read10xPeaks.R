#' Load ATAC-seq data from a 10X Genomics experiment
#' 
#' Creates a \linkS4class{SingleCellExperiment} from the CellRanger output directories for 10X Genomics ATAC-seq data.
#' 
#' @param sample A string containing the path to a single directory corresponding to the output of a CellRanger run.
#' Alternatively, each string may contain a path to a HDF5 file in the sparse matrix format peakrated by 10X.
#' @param type String specifying the type of 10X format to read data from.
#' 
#' @return A \linkS4class{SingleCellExperiment} object containing count data for each peak/motif (row) and cell (column).
#' \itemize{
#' \item Row metadata will either be a \code{\link{GRanges}} object if the features are peaks,
#' or a \linkS4class{DataFrame} of transcription factor motif identifiers and names.
#' \item Column metadata will contain the \code{"Barcode"} field, the cell barcode sequence for each cell library. 
#' \item Rows are named with the peak identifier, while columns are named with the cell barcode.
#' \item The assays will contain a single \code{"counts"} matrix, containing UMI counts for each peak/motif in each cell.
#' Note that the matrix representation will depend on the format of the \code{sample}, see Details.
#' }
#' 
#' @details
#' If \code{type="auto"}, the format of the input file is automatically detected for each \code{samples} based on whether it ends with \code{".h5"}.
#' If so, \code{type} is set to \code{"HDF5"}; otherwise it is set to \code{"sparse"}.
#' \itemize{
#' \item If \code{type="sparse"}, count data are loaded as a \linkS4class{dgCMatrix} object.
#' This is a conventional column-sparse compressed matrix format produced by the CellRanger pipeline,
#' consisting of a (possibly Gzipped) MatrixMarket text file (\code{"matrix.mtx"})
#' with additional tab-delimited files for barcodes (\code{"barcodes.tsv"})
#' and peak/motif annotation (\code{"peaks.bed"} or \code{"motifs.tsv"}).
#' \item If \code{type="prefix"}, count data are also loaded as a \linkS4class{dgCMatrix} object.
#' This assumes the same three-file structure for each sample as described for \code{type="sparse"},
#' but each sample is defined here by a prefix in the file names rather than by being a separate directory.
#' For example, if the \code{samples} entry is \code{"xyx_"},
#' the files are expected to be \code{"xyz_matrix.mtx"}, \code{"xyz_barcodes.tsv"}, etc.
#' \item If \code{type="HDF5"}, count data are assumed to follow the 10X sparse HDF5 format for large data sets.
#' It is loaded as a \linkS4class{TENxMatrix} object, which is a stub object that refers back to the file in \code{samples}.
#' }
#'
#' When \code{type="sparse"} or \code{"prefix"} and \code{compressed=NULL},
#' the function will automatically search for both the unzipped and Gzipped versions of the files.
#' This assumes that the compressed files have an additional \code{".gz"} suffix.
#' We can restrict to only compressed or uncompressed files by setting \code{compressed=TRUE} or \code{FALSE}, respectively.
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{read10xCounts}}, to create the equivalent object for scRNA-seq data.
#' 
#' @examples
#' # Mocking up some 10X genomics output.
#' example(write10xCounts)
#' 
#' # Reading it in.
#' sce10x <- read10xCounts(tmpdir)
#' 
#' # Column names are dropped with multiple 'samples'.
#' sce10x2 <- read10xCounts(c(tmpdir, tmpdir))
#' 
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment 
#' @importFrom SummarizedExperiment rowData<- rowRanges<-
read10xPeaks <- function(sample, type=c("auto", "sparse", "HDF5","prefix"), compressed=NULL) {
    type <- match.arg(type)
    cur.type <- .type_chooser(sample, type)

    if (cur.type=="sparse") {
        info <- .read_peaks_from_sparse(sample, compressed=compressed)
    } else if (cur.type=="prefix") {
        info <- .read_peaks_from_sparse(sample, is.prefix=TRUE, compressed=compressed)
    } else {
        info <- .read_peaks_from_hdf5(sample)
    }

    cell_info <- DataFrame(Barcode = info$cell.names, row.names=info$cell.names)
    sce <- SingleCellExperiment(list(counts = info$mat), colData = cell_info)

    if (is(info$feat.info, "GRanges")) {
        rowRanges(sce) <- info$feat.info
    } else {
        rowData(sce) <- info$feat.info
    }

    sce
}

#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Matrix readMM
#' @importFrom utils read.delim
.read_peaks_from_sparse <- function(path, is.prefix=FALSE, compressed=NULL) {
    FUN <- if (is.prefix) paste0 else file.path

    bname <- "barcodes.tsv"
    mname <- "matrix.mtx"

    barcode.loc <- FUN(path, bname)
    matrix.loc <- FUN(path, mname)

    barcode.loc <- .check_for_compressed(barcode.loc, compressed)
    matrix.loc <- .check_for_compressed(matrix.loc, compressed)

    fname <- "peaks.bed"
    feat.loc <- FUN(path, fname)
    feat.loc <- .check_for_compressed(feat.loc, compressed, error=FALSE)

    if (!file.exists(feat.loc)) {
        feat.loc <- "motifs.tsv"
        feat.loc <- FUN(path, fname)
        feat.loc <- .check_for_compressed(feat.loc, compressed, error=FALSE)

        if (!file.exists(feat.loc)) {
            stop("cannot find 'peaks.bed' or 'motifs.tsv'")
        }

        feat.info <- read.delim(feat.loc)
        colnames(feat.info) <- c("ID", "Name")
    } else {
        feat.info <- rtracklayer::import(feat.loc)
    }

    list(
        mat=as(readMM(matrix.loc), "dgCMatrix"),
        cell.names=readLines(barcode.loc),
        feat.info=feat.info
    )
}

#' @importFrom rhdf5 h5read
#' @importFrom HDF5Array TENxMatrix
#' @importFrom S4Vectors DataFrame
.read_peaks_from_hdf5 <- function(path, genome=NULL, version) {
    group <- "matrix"

    peak.info <- as.character(h5read(path, paste0(group, "/features/id")))
    ftypes <- as.character(h5read(path, paste0(group, "/features/feature_type")))
    if (all(ftypes=="Peaks")) {
        peak.info <- GenomicRanges::GRanges(peak.info)
    } else {
        peak.info <- DataFrame(row.names=peak.info, Feature=ftypes)
    }

    mat <- TENxMatrix(path, group)
    dimnames(mat) <- NULL # for consistency.

    list(
        mat=mat,
        cell.names=as.character(h5read(path, paste0(group, "/barcodes"))),
        peak.info=peak.info
    )
}
