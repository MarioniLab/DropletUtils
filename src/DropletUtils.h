#ifndef DROPLETUTILS_H
#define DROPLETUTILS_H

#include "Rcpp.h"
#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"

#include <algorithm>
#include <deque>
#include <stdexcept>
#include <functional>
#include <cmath>

// Functions to be called from R.

extern "C" {

SEXP downsample_matrix(SEXP, SEXP);

SEXP calculate_random_dev(SEXP, SEXP);

SEXP calculate_pval(SEXP, SEXP, SEXP, SEXP, SEXP);

}

#endif

