#include "DropletUtils.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/utils/const_column.h"
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
SEXP compute_multinom_internal(SEXP mat, SEXP prop, SEXP alpha) {
    auto M=beachmat::create_matrix<MAT>(mat);
    const size_t NC=M->get_ncol();
    const size_t NR=M->get_nrow();
    Rcpp::NumericVector output(NC);

    Rcpp::NumericVector ambient(prop);
    if (ambient.size()!=NR) {
        throw std::runtime_error("length of ambient vector should be equal to number of columns");
    }
    likelihood_calculator lc(alpha);
    beachmat::const_column<MAT> col_holder(M.get());

    for (size_t c=0; c<NC; ++c) {
        double& cur_out=output[c];
        col_holder.fill(c);
        auto val=col_holder.get_values();

        if (col_holder.is_sparse()) {
            size_t num=col_holder.get_n();
            auto dex=col_holder.get_indices();
            for (size_t i=0; i<num; ++i, ++val, ++dex) {
                cur_out+=lc(*val, ambient[*dex]);
            }

        } else {
            for (auto aIt=ambient.begin(); aIt!=ambient.end(); ++aIt, ++val) {
                if (*val) { cur_out+=lc(*val, *aIt); }
            }
        }
    }

    return output;
}

SEXP compute_multinom(SEXP mat, SEXP prop, SEXP alpha) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(mat);
    if (rtype==INTSXP) {
        return compute_multinom_internal<Rcpp::IntegerVector, beachmat::integer_matrix>(mat, prop, alpha);
    } else {
        return compute_multinom_internal<Rcpp::NumericVector, beachmat::numeric_matrix>(mat, prop, alpha);
    }
    END_RCPP    
}

