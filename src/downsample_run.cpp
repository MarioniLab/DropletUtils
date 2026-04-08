#include "Rcpp.h"
#include "scuttle/downsample.h"
#include <stdexcept>

template<typename Input_>
double sum(Input_ start, Input_ end) {
    double total = 0;
    for (auto it = start; it != end; ++it) {
        const auto r = *it;
        if (r < 0) {
            throw std::runtime_error("number of reads should be non-negative");
        }
        total += r;
    }
    if (scuttle::too_large_for_integer_precision(total)) {
        throw std::runtime_error("total number of reads is too large for integer precision");
    } 
    return total;
}

//[[Rcpp::export]]
Rcpp::IntegerVector downsample_run(Rcpp::IntegerVector reads, double prop) {
    double total = sum(reads.begin(), reads.end());
    if (prop < 0 || prop > 1) {
        throw std::runtime_error("'prop' should be in [0, 1]");
    }
    double required = total * prop;

    Rcpp::IntegerVector output(reads.size());
    scuttle::downsample(reads.begin(), reads.end(), output.begin(), total, required);
    return output;
}

//[[Rcpp::export]]
Rcpp::IntegerVector downsample_run_per_cell(Rcpp::IntegerVector cells, Rcpp::IntegerVector reads, Rcpp::NumericVector prop) {
    if (cells.size() != prop.size()) {
        throw std::runtime_error("'cells' and 'prop' should be of the same length");
    }

    // Setting up the output.
    Rcpp::IntegerVector output(reads.size());
    auto oIt = output.begin();
    auto rIt = reads.begin();
    auto pIt = prop.begin();

    // Iterating across the molecule cells and downsampling.
    for (const auto& cell : cells) {
        double total = sum(rIt, rIt + cell);
        const auto prop = *pIt;
        if (prop < 0 || prop > 1) {
            throw std::runtime_error("'prop' should be in [0, 1]");
        }
        double required = total * prop;

        scuttle::downsample(rIt, rIt + cell, oIt, total, required);
        ++pIt;
        rIt += cell;
        oIt += cell;
    }

    return output;
}
