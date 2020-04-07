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

    const double mean_prop=std::accumulate(prop.begin(), prop.end(), 0)/NR; 
    const int upto=std::min(NR, 3);

    Rcpp::IntegerVector output_best(NC), output_second(NC);
    Rcpp::NumericVector output_fc(NC), output_fc2(NC);
    V tmp(NR);
    std::vector<std::pair<double, int> > collected(NR);

    for (int i=0; i<NC; ++i) {
        M->get_col(i, tmp.begin());
        for (int j=0; j<NR; ++j) {
            collected[j].first=tmp[j]/prop[j];
            collected[j].second=j;
        }

        // Estimating the scaling effect.
        const int half=NR/2;
        std::partial_sort(collected.begin(), collected.begin() + half + 1, 
            collected.end(), std::greater<std::pair<double, int> >());

        double scaling;
        if (NR <= 3) {
            scaling=collected[NR-1].first;
        } else if (NR == 4) {
            scaling=collected[2].first;
        } else if ((NR % 2)==1) {
            scaling=collected[half].first;
        } else {
            scaling=(collected[half-1].first + collected[half].first)/2;
        }

        /* First, subtracting the ambient solution. Then adding back a constant
         * to be used as the pseudo-count. The product `scaling * mean_prop`
         * represents the average ambient contamination for each HTO, so we just
         * use it as a pseudo-count (scaled by the `pseudo_count` value). Nice 
         * thing about this is that it adapts to the sequencing depth, under 
         * some assumptions; see documentation for more details.
         *
         * We do require, at least, that the pseudo count be at least
         * `pseudo_count`, just to avoid problems at `scaling=0`!
         */
        const double PSEUDO=std::max(1.0, scaling * mean_prop) * pseudo_count; 

        for (auto& x : collected) {
            x.first = tmp[x.second] - scaling * prop[x.second];
            x.first = std::max(x.first, 0.0) + PSEUDO;
        }

        // Estimating the actual log-fold changes now. 
        std::partial_sort(collected.begin(), collected.begin()+upto,
            collected.end(), std::greater<std::pair<double, int> >());

        if (upto < 1) {
            output_best[i]=NA_INTEGER;
        } else {
            output_best[i]=collected[0].second;
        }

        if (upto < 2) {
            output_fc[i]=R_NaReal;
        } else {
            output_fc[i]=collected[0].first/collected[1].first;
        }

        if (upto < 3) {
            output_second[i]=NA_INTEGER;
            output_fc2[i]=R_NaReal;
        } else {
            output_second[i]=collected[1].second;
            output_fc2[i]=collected[1].first/collected[2].first;
        }
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
