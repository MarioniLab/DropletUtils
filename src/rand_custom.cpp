#include "rand_custom.h"

#include <sstream>
#include <stdexcept>

void check_pcg_vectors(const Rcpp::NumericVector& seeds, const Rcpp::IntegerVector& streams, size_t N, const char* msg) {
    // Use of references to Rcpp classes is deliberate, to ensure that seeds are NumericVectors 
    // and not coerced to integer at some point (which would give NAs or UB for [0, 2^32) integers
    // outside R's maximum range.

    if (seeds.size()!=N) {
        std::stringstream err;
        err << "number of " << msg << " and seeds should be the same";
        throw std::runtime_error(err.str());
    }

    if (streams.size()!=N) {
        std::stringstream err;
        err << "number of " << msg << " and streams should be the same";
        throw std::runtime_error(err.str());
    }

    return;
}

pcg32 create_pcg32(const Rcpp::NumericVector& seeds, const Rcpp::IntegerVector& streams, size_t i) {
    const double curseed=seeds[i];
    if (curseed < 0 || curseed >= 4294967296.0) {
        throw std::runtime_error("seed for PCG32 must lie in [0, 2^32)");
    }

    const double curstream=streams[i];
    if (curstream < 0) { // no need to check upper bound for 32-bit signed ints.
        throw std::runtime_error("stream for PCG32 must be non-negative");
    }

    return pcg32(static_cast<uint32_t>(curseed), curstream);
}
