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
    std::vector<std::pair<double, int> > corrected;

    hash_ambient_adjuster(const Rcpp::NumericVector& p, int pseudo, int nexp) : 
        NR(p.size()), half(NR/2), 
        prop(p), mean_prop(std::accumulate(prop.begin(), prop.end(), 0.0)/NR),
        pseudo_count(pseudo), n_expected(nexp),
        corrected(NR) 
    {
        for (int j=0; j<NR; ++j) {
            if (prop[j] <= 0 || !R_FINITE(prop[j])) {
                throw std::runtime_error("'prop' should only contain positive values");
            }
        }
    }

    double correct(const double* ptr) {
        if (NR==0) {
            return 0;
        }

        for (int j=0; j<NR; ++j) {
            corrected[j].first=ptr[j]/prop[j];
            corrected[j].second=j;
        }

        double scaling=0;

        /* We use either the smallest ratio or the one that is twice the
         * number of expected genuine HTOs (to allow for doublets). Note,
         * no -1 on the latter, we want the first one _past_ the doublet.
         */
        const int fallback = std::min(NR - 1, 2 * n_expected); 

        if (half > fallback) {
            // Computing the actual median.
            std::partial_sort(corrected.begin(), corrected.begin() + half + 1, 
                corrected.end(), std::greater<std::pair<double, int> >());

            if ((NR % 2) == 0) {
                scaling = (corrected[half-1].first + corrected[half].first)/2;
            } else {
                scaling = corrected[half].first;
            }
        } else {
            std::partial_sort(corrected.begin(), corrected.begin() + fallback + 1, 
                corrected.end(), std::greater<std::pair<double, int> >());
            scaling = corrected[fallback].first;
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

        for (auto& x : corrected) {
            const double expected = scaling * prop[x.second];
            const double& observed = ptr[x.second];

            /* 'observed - expected' are the actual counts for this HTO in
             * this cell, then we add back 'PSEUDO' so that every HTO has 
             * exactly the average amount of ambient contamination.
             */
            x.first = std::max(0.0, observed - expected) + PSEUDO;
        }

        // Finally, resorting so that the corrected values are in decreasing order.
        const int upto=std::min(NR, n_expected + 1);
        std::partial_sort(corrected.begin(), corrected.begin()+upto,
            corrected.end(), std::greater<std::pair<double, int> >());

        return PSEUDO;
    }

    template<class IT>
    void fill_best(IT start) {
        if (NR < n_expected) {
            std::fill(start, start + n_expected, NA_INTEGER);
        } else {
            auto ciIt = start;
            for (int e = 0; e < n_expected; ++e, ++ciIt) {
                *ciIt = corrected[e].second;
            }
            std::sort(start, start + n_expected);
        }
    }
};

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

    for (int i=0; i<NC; ++i) {
        auto ptr = M->get_col(i, tmp.data());
        auto PSEUDO = hash.correct(ptr);

        auto COL = output_best.column(i);
        hash.fill_best(COL.begin());

        auto& corrected=hash.corrected;
        if (NR < n_expected + 1) {
            output_fc[i] = R_NaReal;
        } else {
            output_fc[i] = corrected[n_expected - 1].first / corrected[n_expected].first;
        }

        if (NR <= n_expected * 2) {
            output_second[i]=NA_INTEGER;
            output_fc2[i]=R_NaReal;
        } else {
            // We use PSEUDO as this is the level of ambient contamination
            // after all adjustments have been applied.
            output_second[i]=corrected[n_expected].second;
            output_fc2[i]=corrected[n_expected].first/PSEUDO;
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
Rcpp::List hashed_constant(Rcpp::RObject mat, Rcpp::NumericVector prop, double pseudo_count, int n_expected) {
    auto M = beachmat::read_lin_block(mat);
    const int NR = M->get_nrow();
    const int NC = M->get_ncol();
    if (prop.size() != NR) {
        throw std::runtime_error("'length(prop)' should be the same as 'nrow(mat)'");
    }

    Rcpp::IntegerMatrix output_best(n_expected, NC); 
    Rcpp::IntegerVector output_second(NC);
    Rcpp::NumericVector output_fc(NC), output_num(NC), output_amb(NC);
    std::vector<double> tmp(NR);

    hash_ambient_adjuster hash(prop, pseudo_count, n_expected);

    for (int i=0; i<NC; ++i) {
        auto ptr = M->get_col(i, tmp.data());
        auto PSEUDO = hash.correct(ptr);

        auto COL = output_best.column(i);
        hash.fill_best(COL.begin());

        auto& corrected=hash.corrected;
        if (NR < n_expected + 1) {
            output_fc[i] = R_NaReal;
        } else {
            output_fc[i] = corrected[n_expected - 1].first / corrected[n_expected].first;
        }

        if (NR <= n_expected * 2) {
            output_second[i]=NA_INTEGER;
            output_num[i]=R_NaReal;
        } else {
            // We use PSEUDO as this is the level of ambient contamination
            // after all adjustments have been applied.
            output_second[i]=corrected[n_expected].second;
            output_num[i]=corrected[n_expected].first;
        }
        output_amb[i]=PSEUDO;
    }

    return Rcpp::List::create(
        Rcpp::Named("Best")=output_best, 
        Rcpp::Named("Second")=output_second,
        Rcpp::Named("FC")=output_fc,
        Rcpp::Named("Numerator")=output_num,
        Rcpp::Named("Ambient")=output_amb
    );
}


