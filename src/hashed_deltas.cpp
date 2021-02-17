#include "Rcpp.h"
#include "beachmat3/beachmat.h"

#include <stdexcept>
#include <algorithm>
#include <vector>

struct hash_ambient_adjuster {
    const int NR, half;
    const Rcpp::NumericVector& prop;
    const double mean_prop;
    const double pseudo_count;
    const int n_expected;
    std::vector<std::pair<double, int> > collected;

    hash_ambient_adjuster(const Rcpp::NumericVector& p, int pseudo, int nexp) : 
        NR(p.size()), half(NR/2), 
        prop(p), mean_prop(std::accumulate(prop.begin(), prop.end(), 0.0)/NR),
        pseudo_count(pseudo), n_expected(nexp),
        collected(NR) 
    {
        for (int j=0; j<NR; ++j) {
            if (prop[j] <= 0 || !R_FINITE(prop[j])) {
                throw std::runtime_error("'prop' should only contain positive values");
            }
        }
    }

    double correct(const double* ptr, bool full=false) {
        if (NR==0) {
            return 0;
        }

        for (int j=0; j<NR; ++j) {
            collected[j].first=ptr[j]/prop[j];
            collected[j].second=j;
        }

        double scaling=0;

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

        /* 'scaling * mean_prop' represents the average ambient contamination
         * for this cell, where the average is taken across all HTOs. We will
         * be adjusting the counts so that the ambient contamination for all
         * HTOs is equal to this average. In effect, the ambient contamination
         * serves as a pseudo-count for the later log-fold change calculations,
         * which has some nice properties as it naturally adjusts to the 
         * sequencing depth - see the documentation for some more details.
         *
         * We do require, at least, that the pseudo count be at least
         * `pseudo_count`, just to avoid problems at `scaling=0`!
         */
        const double PSEUDO=std::max(pseudo_count, scaling * mean_prop); 

        for (auto& x : collected) {
            const double expected = scaling * prop[x.second];
            const double& observed = ptr[x.second];

            /* 'observed - expected' are the actual counts for this HTO in
             * this cell, then we add back 'PSEUDO' so that every HTO has 
             * exactly the average amount of ambient contamination.
             */
            x.first = std::max(0.0, observed - expected) + PSEUDO;
        }

        return PSEUDO;
    }
};

// [[Rcpp::export(rng=false)]]
Rcpp::List hashed_subtract(Rcpp::RObject mat, Rcpp::NumericVector prop, double pseudo_count, int n_expected) {
    auto M = beachmat::read_lin_block(mat);
    const int NR = M->get_nrow();
    const int NC = M->get_ncol();
    if (prop.size() != NR) {
        throw std::runtime_error("'length(prop)' should be the same as 'nrow(mat)'");
    }

    Rcpp::NumericVector outp(NC);
    Rcpp::NumericMatrix output(NR, NC); 
    std::vector<double> tmp(NR);

    hash_ambient_adjuster hash(prop, pseudo_count, n_expected);

    for (int i=0; i<NC; ++i) {
        auto ptr = M->get_col(i, tmp.data());
        outp[i] = hash.correct(ptr);

        auto& collected=hash.collected;
        auto outcol = output.column(i);
        for (const auto& pv : collected) {
            outcol[pv.second] = pv.first;
        }
    }

    return Rcpp::List::create(outp, output);
}

// [[Rcpp::export(rng=false)]]
Rcpp::List hashed_deltas(Rcpp::RObject mat, Rcpp::NumericVector prop, double pseudo_count, int n_expected) {
    auto M = beachmat::read_lin_block(mat);
    const int NR = M->get_nrow();
    const int NC = M->get_ncol();
    if (prop.size() != NR) {
        throw std::runtime_error("'length(prop)' should be the same as 'nrow(mat)'");
    }

    Rcpp::IntegerMatrix output_best(n_expected, NC); 
    Rcpp::IntegerVector output_second(NC);
    Rcpp::NumericVector output_fc(NC), output_fc2(NC);
    std::vector<double> tmp(NR);

    hash_ambient_adjuster hash(prop, pseudo_count, n_expected);
    auto& collected=hash.collected;

    for (int i=0; i<NC; ++i) {
        auto ptr = M->get_col(i, tmp.data());
        auto PSEUDO=hash.correct(ptr);

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

        if (upto <= n_expected * 2) {
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


