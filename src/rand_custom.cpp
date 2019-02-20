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

pcg32 create_pcg32(SEXP seeds, int stream) {
    return pcg32(dqrng::convert_seed<uint64_t>(seeds), stream);
}
