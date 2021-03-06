# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

downsample_run <- function(reads, prop) {
    .Call('_DropletUtils_downsample_run', PACKAGE = 'DropletUtils', reads, prop)
}

downsample_run_per_cell <- function(cells, reads, prop) {
    .Call('_DropletUtils_downsample_run_per_cell', PACKAGE = 'DropletUtils', cells, reads, prop)
}

encode_sequences <- function(Seqs) {
    .Call('_DropletUtils_encode_sequences', PACKAGE = 'DropletUtils', Seqs)
}

find_chimeric <- function(cells, umis, reads, minfrac, diagnostics) {
    .Call('_DropletUtils_find_chimeric', PACKAGE = 'DropletUtils', cells, umis, reads, minfrac, diagnostics)
}

find_swapped <- function(cells, genes, umis, reads, minfrac, diagnostics) {
    .Call('_DropletUtils_find_swapped', PACKAGE = 'DropletUtils', cells, genes, umis, reads, minfrac, diagnostics)
}

get_cell_barcodes <- function(fname, dname, barcodelen) {
    .Call('_DropletUtils_get_cell_barcodes', PACKAGE = 'DropletUtils', fname, dname, barcodelen)
}

group_cells <- function(cells, gems) {
    .Call('_DropletUtils_group_cells', PACKAGE = 'DropletUtils', cells, gems)
}

hashed_deltas <- function(mat, prop, pseudo_count, n_expected) {
    .Call('_DropletUtils_hashed_deltas', PACKAGE = 'DropletUtils', mat, prop, pseudo_count, n_expected)
}

hashed_constant <- function(mat, prop, pseudo_count, n_expected) {
    .Call('_DropletUtils_hashed_constant', PACKAGE = 'DropletUtils', mat, prop, pseudo_count, n_expected)
}

montecarlo_pval <- function(totalval, totallen, prob, ambient, iterations, alpha, seeds, streams) {
    .Call('_DropletUtils_montecarlo_pval', PACKAGE = 'DropletUtils', totalval, totallen, prob, ambient, iterations, alpha, seeds, streams)
}

