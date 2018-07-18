#include "DropletUtils.h"

template <typename V, class MAT>
SEXP compute_multinom_internal(MAT M, SEXP prop, SEXP alpha) {
    const size_t NC=M->get_ncol();
    const size_t NR=M->get_nrow();
    V target(NR);
    Rcpp::NumericVector output(NC);

    Rcpp::NumericVector ambient(prop);
    if (ambient.size()!=NR) {
        throw std::runtime_error("length of ambient vector should be equal to number of columns");
    }

    double Alpha=check_numeric_scalar(alpha, "alpha");
    const bool use_alpha=R_FINITE(Alpha);
    if (use_alpha && Alpha <= 0) {
        throw std::runtime_error("alpha must be positive if specified");
    }

    for (size_t c=0; c<NC; ++c) {
        auto info=M->get_const_col_indexed(c, target.begin());
        size_t num=std::get<0>(info);
        auto dex=std::get<1>(info);
        auto val=std::get<2>(info);
        
        double& cur_out=output[c];
        for (size_t i=0; i<num; ++i, ++val, ++dex) {
            if (use_alpha) {
                cur_out += std::lgamma(*val + ambient[*dex] * Alpha) - std::lgamma(*val + 1) - std::lgamma(ambient[*dex] * Alpha);
            } else {
                cur_out += (*val) * std::log(ambient[*dex]) - std::lgamma(*val + 1);
            }
        }
    }

    return output;
}

SEXP compute_multinom(SEXP mat, SEXP prop, SEXP alpha) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(mat);
    if (rtype==INTSXP) {
        auto ptr=beachmat::create_integer_matrix(mat);
        return compute_multinom_internal<Rcpp::IntegerVector>(ptr.get(), prop, alpha);
    } else {
        auto ptr=beachmat::create_numeric_matrix(mat);
        return compute_multinom_internal<Rcpp::NumericVector>(ptr.get(), prop, alpha);
    }
    END_RCPP    
}

