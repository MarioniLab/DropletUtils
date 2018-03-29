#include "DropletUtils.h"

template<class V, class M>
SEXP load_tenx_to_hdf5_internal(SEXP fhandle, SEXP chunksize, SEXP x_type, SEXP nz, M omat) {
    const int chunks=check_integer_scalar(chunksize, "chunk size");
    const size_t NZ=check_integer_scalar(nz, "number of non-zero entries");
    
    Rcpp::Environment pkgenv=Rcpp::Environment::namespace_env("DropletUtils");
    Rcpp::Function scanfun=pkgenv[".scan_mm_file"];
    Rcpp::LogicalVector reorder(1, 1);

    Rcpp::IntegerVector holding(chunks);
    size_t counter=0;

    while (counter < NZ) {
        Rcpp::List incoming=scanfun(fhandle, chunksize, x_type, reorder);
        if (incoming.size()!=3) {
            throw std::runtime_error("scan() should return a list of three integer vectors");
        }
        Rcpp::IntegerVector I=incoming[0], J=incoming[1];
        V X=incoming[2];
        if (I.size()!=J.size() || I.size()!=X.size()) {
            throw std::runtime_error("integer vectors should be of the same length");
        }

        // Iterating through; pulling the column out, editing it, and then saving the column again.
        auto IIt=I.begin();
        auto JIt=J.begin();
        auto XIt=X.begin();

        while (JIt!=J.end()) {
            int curcol=*JIt;
            size_t counter=0;
            auto hIt=holding.begin();

            while (JIt!=J.end() && curcol==*JIt) {
                (*hIt)=*IIt - 1; // zero indexing.
                ++counter;
                ++JIt;
                ++IIt;
                ++hIt;
            }

            omat->set_col_indexed(curcol-1, counter, holding.begin(), XIt);
            XIt+=counter;
        }

        counter += chunks;
    }

    return omat->yield();
}

SEXP load_tenx_to_hdf5(SEXP fhandle, SEXP chunksize, SEXP x_type, SEXP nr, SEXP nc, SEXP nz) {
    BEGIN_RCPP
    const size_t NR=check_integer_scalar(nr, "number of rows");
    const size_t NC=check_integer_scalar(nc, "number of columns");

    if (Rcpp::RObject(x_type).sexp_type()==INTSXP) {
        auto omat=beachmat::create_integer_output(NR, NC, beachmat::HDF5_PARAM);
        return load_tenx_to_hdf5_internal<Rcpp::IntegerVector>(fhandle, chunksize, x_type, nz, omat.get());
    } else { 
        auto omat=beachmat::create_numeric_output(NR, NC, beachmat::HDF5_PARAM);
        return load_tenx_to_hdf5_internal<Rcpp::NumericVector>(fhandle, chunksize, x_type, nz, omat.get());
    }
    END_RCPP
}
