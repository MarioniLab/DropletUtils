#include "Rcpp.h"

#include "boost/random.hpp"
#include "rand_custom.h"
#include "utils.h"

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <cmath>

//[[Rcpp::export(rng=false)]]
Rcpp::IntegerVector montecarlo_pval (Rcpp::IntegerVector totalval, Rcpp::IntegerVector totallen, 
    Rcpp::NumericVector prob, Rcpp::NumericVector ambient, int iterations, double alpha, 
    Rcpp::List seeds, Rcpp::IntegerVector streams) 
{
    if (totalval.size()!=totallen.size()) {
        throw std::runtime_error("length of run value/length vectors are not the same");
    }
    const int nvalues=std::accumulate(totallen.begin(), totallen.end(), 0);
    if (nvalues!=prob.size()) {
        throw std::runtime_error("sum of run lengths does not equal length of 'P' vector");
    }

    if (iterations < 0) {
        throw std::runtime_error("number of iterations should be a non-negative integer");
    }
    check_pcg_vectors(seeds, streams, iterations, "iterations");

    const bool use_alpha=R_FINITE(alpha);
    if (use_alpha && alpha <= 0) {
        throw std::runtime_error("alpha must be positive if specified");
    }

    // Setting up temporary objects.
    Rcpp::IntegerVector above(nvalues);
    const size_t ngenes=ambient.size();
    if (ngenes==0) {
        return above;
    }

    std::vector<int> tracker(ngenes);
    std::vector<double> probs(ngenes), logprob;
    if (use_alpha) {
        probs.resize(ngenes);
    } else {
        logprob.reserve(ngenes);
        for (auto a : ambient) {
            logprob.push_back(std::log(a));
        }
        std::copy(ambient.begin(), ambient.end(), probs.begin());
    }

    // Looping across iterations, using a new probability vector per iteration.
    for (int it=0; it<iterations; ++it) {
        auto generator=create_pcg32(seeds[it], streams[it]);

        if (use_alpha) {
            typedef boost::random::gamma_distribution<double> distr_t;
            typedef typename distr_t::param_type param_t;

            // Do NOT cache gamma distribution across iterations, as this introduces possible dependencies.
            distr_t cpp_gamma;
            for (size_t ldx=0; ldx<ngenes; ++ldx) {
                probs[ldx]=cpp_gamma(generator, param_t(ambient[ldx] * alpha, 1));
            }
        }
        boost::random::discrete_distribution<> sampler(probs.begin(), probs.end());

        int curtotal=0;
        std::fill(tracker.begin(), tracker.end(), 0);
        double curp=0;

        auto abIt=above.begin();
        auto tvIt=totalval.begin();
        auto tlIt=totallen.begin();
        auto pIt=prob.begin();

        while (tvIt!=totalval.end()) {
            const auto& curlen=*tlIt;
            const auto& curval=*tvIt;

            // Sampling more points to reach the current total count.
            while (curtotal<curval) {
                auto chosen=sampler(generator);
                auto& curnum=tracker[chosen];
                if (use_alpha) { 
                    // Difference of Gamma's in data-dependent component of Dirichlet-multinomial.
                    curp += std::log(ambient[chosen] * alpha + curnum); 
                } else {
                    curp += logprob[chosen];
                }
                curp -= std::log(++curnum); // corresponds to division by factorial on observed count.
                ++curtotal;
            }

            // Figuring out where it lies in the probability vector.
            size_t higher=std::lower_bound(pIt, pIt+curlen, curp) - pIt;
            if (higher<curlen) {
                ++(*(abIt+higher));
            }

            ++tlIt;
            ++tvIt;
            pIt+=curlen;
            abIt+=curlen;
        }
    }

    // Calculating the number of simulations above each value.
    auto abIt=above.begin();
    for (const auto& curlen : totallen) {
        for (int i=1; i<curlen; ++i) {
            const auto& prev=*abIt;
            ++abIt;
            (*abIt)+=prev;
        }
        ++abIt;
    }

    return above;
}
