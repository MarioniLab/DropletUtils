#include "rand_custom.h"

#include <sstream>
#include <stdexcept>

void check_pcg_vectors(const Rcpp::IntegerVector& seeds, const Rcpp::IntegerVector& streams, size_t N, const char* msg) {
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

pcg32 create_pcg32(const Rcpp::IntegerVector& seeds, const Rcpp::IntegerVector& streams, size_t i) {
    const int curseed=seeds[i];
    if (curseed < 0) {
        throw std::runtime_error("seed for PCG32 must be non-negative");
    }

    const int curstream=streams[i];
    if (curstream < 0) { // no need to check upper bound for 32-bit signed ints.
        throw std::runtime_error("stream for PCG32 must be non-negative");
    }

    return pcg32(curseed, curstream);
}
