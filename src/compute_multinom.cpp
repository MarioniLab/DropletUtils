#include "DropletUtils.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "utils.h"

#include <stdexcept>
#include <cmath>
#include <tuple>

class likelihood_calculator {
public: 
    likelihood_calculator(SEXP alpha) : Alpha(check_numeric_scalar(alpha, "alpha")), use_alpha(use_alpha=R_FINITE(Alpha)) {
        if (use_alpha && Alpha <= 0) {
            throw std::runtime_error("alpha must be positive if specified");
        }
    }

    double operator()(double count, double prop) const {
        if (use_alpha) {
            return std::lgamma(count + prop * Alpha) - std::lgamma(count + 1) - std::lgamma(prop * Alpha);
        } else {
            return count * std::log(prop) - std::lgamma(count + 1);
        }
    }
private:
    double Alpha;
    bool use_alpha;
};


template <typename V, class MAT>
SEXP compute_multinom_internal(MAT M, SEXP prop, SEXP alpha) {
    const size_t NC=M->get_ncol();
    const size_t NR=M->get_nrow();
    Rcpp::NumericVector output(NC);

    Rcpp::NumericVector ambient(prop);
    if (ambient.size()!=NR) {
        throw std::runtime_error("length of ambient vector should be equal to number of columns");
    }
    likelihood_calculator lc(alpha);

    const bool is_sparse=M->col_raw_type()=="sparse";
    const bool is_dense=M->col_raw_type()=="dense";
    V target(NR);
    auto raws=M->set_up_raw();

    for (size_t c=0; c<NC; ++c) {
        double& cur_out=output[c];

        if (is_sparse) {
            M->get_col_raw(c, raws);
            size_t num=raws.get_n();
            auto dex=raws.get_structure_start();
            auto val=raws.get_values_start();
            for (size_t i=0; i<num; ++i, ++val, ++dex) {
                cur_out+=lc(*val, ambient[*dex]);
            }

        } else {
            typename V::iterator it;
            if (is_dense) {
                M->get_col_raw(c, raws);
                it=raws.get_values_start();
            } else {
                it=target.begin();
                M->get_col(c, it);
            }

            auto aIt=ambient.begin();
            while (aIt!=ambient.end()) {
                if (*it) { cur_out+=lc(*it, *aIt); }
                ++it;
                ++aIt;
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

