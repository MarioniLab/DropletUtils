#include "Rcpp.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"

#include <stdexcept>
#include <algorithm>
#include <vector>

template <typename V, class MAT>
Rcpp::List hashed_deltas_internal(Rcpp::RObject mat, Rcpp::NumericVector prop, double pseudo_count) {
    auto M=beachmat::create_matrix<MAT>(mat);
    const int NR=M->get_nrow();
    const int NC=M->get_ncol();
    if (prop.size()!=NR) {
        throw std::runtime_error("'length(prop)' should be the same as 'nrow(mat)'");
    }
    if (NR<3) {
        throw std::runtime_error("'nrow(mat)' should be at least 3 for doublet detection");
    }

    Rcpp::IntegerVector output_best(NC), output_second(NC);
    Rcpp::NumericVector output_fc(NC), output_fc2(NC);
    V tmp(NR);
    std::vector<std::pair<double, int> > collected(NR);

    for (int i=0; i<NC; ++i) {
        M->get_col(i, tmp.begin());
        for (int j=0; j<NR; ++j) {
            collected[j].first=tmp[j]/prop[j] + pseudo_count;
            collected[j].second=j;
        }

        std::partial_sort(collected.begin(), collected.begin()+3,
            collected.end(), std::greater<std::pair<double, int> >());
        output_best[i]=collected[0].second;
        output_second[i]=collected[1].second;
        output_fc[i]=collected[0].first/collected[1].first;
        output_fc2[i]=collected[1].first/collected[2].first;
    }

    return Rcpp::List::create(
        Rcpp::Named("Best")=output_best, 
        Rcpp::Named("Second")=output_second,
        Rcpp::Named("FC")=output_fc,
        Rcpp::Named("FC2")=output_fc2
    );
}

// [[Rcpp::export(rng=false)]]
Rcpp::List hashed_deltas(Rcpp::RObject mat, Rcpp::NumericVector prop, double pseudo_count) 
{
    int rtype=beachmat::find_sexp_type(mat);
    if (rtype==INTSXP) {
        return hashed_deltas_internal<Rcpp::IntegerVector, beachmat::integer_matrix>(mat, prop, pseudo_count);
    } else {
        return hashed_deltas_internal<Rcpp::NumericVector, beachmat::numeric_matrix>(mat, prop, pseudo_count);
    }
}
