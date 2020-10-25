#include "Rcpp.h"
#include "beachmat3/beachmat.h"

#include <stdexcept>
#include <algorithm>
#include <vector>

// [[Rcpp::export(rng=false)]]
Rcpp::List hashed_deltas(Rcpp::RObject mat, Rcpp::NumericVector prop, double pseudo_count, int n_expected) {
    auto M = beachmat::read_lin_block(mat);
    const int NR = M->get_nrow();
    const int NC = M->get_ncol();
    if (prop.size() != NR) {
        throw std::runtime_error("'length(prop)' should be the same as 'nrow(mat)'");
    }

    const double mean_prop=std::accumulate(prop.begin(), prop.end(), 0.0)/NR; 

    Rcpp::IntegerMatrix output_best(n_expected, NC); 
    Rcpp::IntegerVector output_second(NC);
    Rcpp::NumericVector output_fc(NC), output_fc2(NC);
    std::vector<double> tmp(NR);
    std::vector<std::pair<double, int> > collected(NR);

    for (int i=0; i<NC; ++i) {
        auto ptr = M->get_col(i, tmp.data());
        for (int j=0; j<NR; ++j) {
            collected[j].first=ptr[j]/prop[j];
            collected[j].second=j;
        }

        double scaling=0;
        if (NR) {
            const int half = NR/2;

            /* We use either the smallest ratio or the one that is twice the
             * number of expected genuine HTOs (to allow for doublets). Note,
             * no -1 on the latter, we want the first one _past_ the doublet.
             */
            const int fallback = std::min(NR - 1, 2 * n_expected); 

            if (half > fallback) {
                // Computing the actual median.
                std::partial_sort(collected.begin(), collected.begin() + half + 1, 
                    collected.end(), std::greater<std::pair<double, int> >());

                if ((NR % 2) == 0) {
                    scaling = (collected[half-1].first + collected[half].first)/2;
                } else {
                    scaling = collected[half].first;
                }
            } else {
                std::partial_sort(collected.begin(), collected.begin() + fallback + 1, 
                    collected.end(), std::greater<std::pair<double, int> >());
                scaling = collected[fallback].first;
            }
        }

        /* First, subtracting the ambient solution. Then adding back a constant
         * to be used as the pseudo-count. The product `scaling * mean_prop`
         * represents the average ambient contamination for each HTO, so we
         * just use it as a pseudo-count. Nice thing about this is that it
         * adapts to the sequencing depth, under some assumptions; see
         * documentation for more details.
         *
         * We do require, at least, that the pseudo count be at least
         * `pseudo_count`, just to avoid problems at `scaling=0`!
         */
        const double PSEUDO=std::max(pseudo_count, scaling * mean_prop); 

        for (auto& x : collected) {
            x.first = ptr[x.second] - scaling * prop[x.second];
            x.first = std::max(x.first, 0.0) + PSEUDO;
        }

        const int upto=std::min(NR, n_expected + 2);
        std::partial_sort(collected.begin(), collected.begin()+upto,
            collected.end(), std::greater<std::pair<double, int> >());

        auto cur_ids = output_best.column(i);
        if (upto < n_expected) {
            std::fill(cur_ids.begin(), cur_ids.end(), NA_INTEGER);
        } else {
            auto ciIt = cur_ids.begin();
            for (int e = 0; e < n_expected; ++e, ++ciIt) {
                *ciIt = collected[e].second;
            }
            std::sort(cur_ids.begin(), cur_ids.end());
        }

        if (upto < n_expected + 1) {
            output_fc[i] = R_NaReal;
        } else {
            output_fc[i] = collected[n_expected - 1].first / collected[n_expected].first;
        }

        if (upto < n_expected + 2) {
            // Technically, I should set this to NR < n_expected * 2, because
            // the scaling factor is not guaranteed to be properly estimated
            // for doublet estimation if that is true. However, it seems
            // unlikely that you'd get a completely complementary set of
            // barcodes in a doublet when doing a choose(n, n/2) barcode
            // combination scheme, so we'll just take the chance.
            output_second[i]=NA_INTEGER;
            output_fc2[i]=R_NaReal;
        } else {
            // We use PSEUDO as this is the level of ambient contamination
            // after all adjustments have been applied.
            output_second[i]=collected[n_expected].second;
            output_fc2[i]=collected[n_expected].first/PSEUDO;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("Best")=output_best, 
        Rcpp::Named("Second")=output_second,
        Rcpp::Named("FC")=output_fc,
        Rcpp::Named("FC2")=output_fc2
    );
}
