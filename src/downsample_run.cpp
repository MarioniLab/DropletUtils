#include "Rcpp.h"
#include "scuttle/downsample_vector.h"
#include <stdexcept>

//[[Rcpp::export]]
Rcpp::IntegerVector downsample_run(Rcpp::IntegerVector reads, double prop) {
    Rcpp::IntegerVector output(reads.size());
    scuttle::downsample_vector(reads.begin(), reads.end(), output.begin(), prop);
    return output;
}

//[[Rcpp::export]]
Rcpp::IntegerVector downsample_run_per_cell(Rcpp::IntegerVector cells, Rcpp::IntegerVector reads, Rcpp::NumericVector prop) {
    if (cells.size()!=prop.size()) {
        throw std::runtime_error("'cells' and 'prop' should be of the same length");
    }

    // Setting up the output.
    Rcpp::IntegerVector output(reads.size());
    auto oIt=output.begin();
    auto rIt=reads.begin();
    auto pIt=prop.begin();

    // Iterating across the molecule cells and downsampling.
    for (const auto& cell : cells) {
        scuttle::downsample_vector(rIt, rIt+cell, oIt, *pIt);
        ++pIt;
        rIt+=cell;
        oIt+=cell;
    }

    return output;
}
