#' @export
#' @importFrom S4Vectors DataFrame
get10xMolInfoStats <- function(sample, barcode.length=NULL) {
    incoming <- read10xMolInfo(sample, barcode.length, get.umi=FALSE)

    out <- .get_cell_ordering(incoming$data$cell, incoming$data$gem_group)
    o <- out$order
    group.id <- out$id
    unique.cells <- out$cell
    unique.gems <- out$gem
    run.length <- out$length # also equal to the number of UMIs.

    num.reads <- by(incoming$data$reads[o], group.id, FUN=sum)
    num.genes <- by(incoming$data$gene[o], group.id, FUN=function(x) { sum(!duplicated(x)) })
    DataFrame(cell=unique.cells, gem_group=unique.gems, num.umis=run.length, 
            num.reads=as.integer(num.reads), num.genes=as.integer(num.genes))
}

.get_cell_ordering <- function(cells, gems) 
# Returns "order(cells, gems)", the cell ID after reordering,
# and the cell barcode/gem group for each cell ID.
# The length refers to the number of entries with the same ID.
{
    out <- .Call(cxx_group_cells, cells, gems)
    names(out) <- c("order", "id", "cell", "gem", "length")
    out$id <- rep(seq_along(out$length), out$length)
    out
}
