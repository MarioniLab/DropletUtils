#include "DropletUtils.h"

SEXP load_tenx_to_hdf5(SEXP fhandle, SEXP chunksize, SEXP nr, SEXP nc, SEXP nz) {
    BEGIN_RCPP

    Rcpp::RObject FHANDLE(fhandle);
    const int chunks=check_integer_scalar(chunksize, "chunk size");
    const size_t NR=check_integer_scalar(nr, "number of rows");
    const size_t NC=check_integer_scalar(nc, "number of columns");
    const size_t NZ=check_integer_scalar(nz, "number of non-zero entries");
    
    Rcpp::Environment pkgenv=Rcpp::Environment::namespace_env("DropletUtils");
    Rcpp::Function scanfun=pkgenv[".scan_mm_file"];
    Rcpp::LogicalVector reorder(1, 1);

    auto omat=beachmat::create_integer_output(NR, NC, beachmat::HDF5_PARAM);
    Rcpp::IntegerVector output(NR);
    size_t counter=0;

    while (counter < NZ) {
        Rcpp::List incoming=scanfun(FHANDLE, chunksize, reorder);
        if (incoming.size()!=3) {
            throw std::runtime_error("scan() should return a list of three integer vectors");
        }
        Rcpp::IntegerVector I=incoming[0], J=incoming[1], X=incoming[2];
        if (I.size()!=J.size() || I.size()!=X.size()) {
            throw std::runtime_error("integer vectors should be of the same length");
        }

        // Iterating through; pulling the column out, editing it, and then saving the column again.
        auto IIt=I.begin();
        auto JIt=J.begin();
        auto XIt=X.begin();
        auto oIt=output.begin()-1; // for zero indexing.

        while (JIt!=J.end()) {
            int curcol=*JIt;
            omat->get_col(curcol-1, output.begin());
            while (JIt!=J.end() && curcol==*JIt) {
                *(oIt + *IIt) = *XIt;
                ++XIt;
                ++IIt;
                ++JIt;
            }
            omat->set_col(curcol-1, output.begin());
        }

        counter += chunks;
    }

    return omat->yield();
    END_RCPP
}
