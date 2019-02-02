#include "rand_custom.h"
#include "convert_seed.h"

#include <sstream>
#include <stdexcept>

void check_pcg_vectors(Rcpp::List seeds, Rcpp::IntegerVector streams, size_t N, const char* msg) {
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

pcg32 create_pcg32(Rcpp::List seeds, Rcpp::IntegerVector streams, size_t i) {
    return pcg32(convert_seed<uint64_t>(seeds[i]), streams[i]);
}
