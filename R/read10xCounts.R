#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment
read10xCounts <- function(samples, col.names=FALSE, type=c("auto", "sparse", "HDF5"), 
    version=c("auto", "2", "3"), genome=NULL) 
# Reads in one or more 10X directories in 'samples', and produces
# a SingleCellExperiment object as the output.
#
# written by Davis McCarthy
# modifications by Aaron Lun
# some time ago.
{
    nsets <- length(samples)
    full_data <- vector("list", nsets)
    gene_info_list <- vector("list", nsets)
    cell_info_list <- vector("list", nsets)

    type <- match.arg(type)
    version <- match.arg(version)

    for (i in seq_len(nsets)) { 
        run <- samples[i]
        cur.type <- .type_chooser(run, type)

        if (cur.type=="sparse") {
            info <- .read_from_sparse(run, version=version)
        } else {
            info <- .read_from_hdf5(run, genome=genome, version=version)
        }

        full_data[[i]] <- info$mat
        gene_info_list[[i]] <- info$gene.info
        cell.names <- info$cell.names
        cell_info_list[[i]] <- DataFrame(Sample = rep(run, length(cell.names)), Barcode = cell.names, row.names=NULL)
    }

    # Checking gene uniqueness. 
    if (nsets > 1 && length(unique(gene_info_list)) != 1L) {
        stop("gene information differs between runs")
    }
    gene_info <- gene_info_list[[1]]
    rownames(gene_info) <- gene_info$ID

    # Forming the full data matrix.
    full_data <- do.call(cbind, full_data)

    # Adding the cell data.
    cell_info <- do.call(rbind, cell_info_list)
    if (col.names) {
        if (nsets == 1L) {
            cnames <- cell_info$Barcode
        } else {
            sid <- rep(seq_along(cell_info_list), vapply(cell_info_list, nrow, 1L))
            cnames <- paste0(sid, "_", cell_info$Barcode)
        }
        colnames(full_data) <- cnames
    }

    SingleCellExperiment(list(counts = full_data), rowData = gene_info, colData = cell_info)
}

#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Matrix readMM
#' @importFrom utils read.delim
.read_from_sparse <- function(path, version) {
    if (version=="auto") {
        version <- if (file.exists(file.path(path, "features.tsv.gz"))) "3" else "2"
    }

    if (version=="3") {
        bname <- "barcodes.tsv.gz"
        gname <- "features.tsv.gz"
        mname <- "matrix.mtx.gz"
    } else {
        bname <- "barcodes.tsv"
        gname <- "genes.tsv"
        mname <- "matrix.mtx"
    }

    barcode.loc <- file.path(path, bname)
    gene.loc <- file.path(path, gname)
    matrix.loc <- file.path(path, mname)

    gene.info <- read.delim(gene.loc, header=FALSE, colClasses="character", stringsAsFactors=FALSE, quote="", comment.char="")
    if (version=="3") {
        colnames(gene.info) <- c("ID", "Symbol", "Type")
    } else {
        colnames(gene.info) <- c("ID", "Symbol")
    }

    list(
        mat=as(readMM(matrix.loc), "dgCMatrix"),
        cell.names=readLines(barcode.loc),
        gene.info=gene.info
    )
}

#' @importFrom rhdf5 h5ls h5read
#' @importFrom utils head
.read_from_hdf5 <- function(path, genome=NULL, version) {
    available <- h5ls(path, recursive=FALSE)
    available <- available[available$otype=="H5I_GROUP",]

    if (version=="auto") {
        version <- if ("matrix" %in% available$name) "3" else "2"
    }

    if (version=="2") {
        group <- genome
        if (is.null(group)) {
            if (nrow(available) > 1L) {
                to.see <- head(available$name, 3)
                if (length(to.see)==3L) {
                    to.see[3] <- "..."
                }
                stop("more than one available group (", paste(to.see, collapse=", "), ")")
            } else if (nrow(available) == 0L) {
                stop("no available groups")
            }
            group <- available$name
        }

        gene.info <- data.frame(
            ID=as.character(h5read(path, paste0(group, "/genes"))),
            Symbol=as.character(h5read(path, paste0(group, "/gene_names"))),
            stringsAsFactors=FALSE
        )
    } else {
        group <- "matrix"
        gene.info <- data.frame(
            ID=as.character(h5read(path, paste0(group, "/features/id"))),
            Symbol=as.character(h5read(path, paste0(group, "/features/name"))),
            Type=as.character(h5read(path, paste0(group, "/features/feature_type"))),
            stringsAsFactors=FALSE
        )
    }

    mat <- TENxMatrix(path, group)
    dimnames(mat) <- NULL # for consistency.
    list(
        mat=mat,
        cell.names=as.character(h5read(path, paste0(group, "/barcodes"))),
        gene.info=gene.info
    )
}
