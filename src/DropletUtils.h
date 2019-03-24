#ifndef DROPLETUTILS_H
#define DROPLETUTILS_H

#include "Rcpp.h"

// Functions to be called from R.

extern "C" {

SEXP downsample_matrix(SEXP, SEXP, SEXP);

SEXP downsample_runs(SEXP, SEXP, SEXP, SEXP);


SEXP compute_multinom(SEXP, SEXP, SEXP);

SEXP montecarlo_pval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


SEXP find_swapped(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


SEXP get_cell_barcodes(SEXP, SEXP, SEXP);

SEXP encode_sequences(SEXP);

SEXP group_cells(SEXP, SEXP);

}

#endif

