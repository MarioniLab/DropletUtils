#include "DropletUtils.h"

SEXP montecarlo_pval (SEXP _totalval, SEXP _totallen, SEXP _P, SEXP _ambient, SEXP _iter) { 
    BEGIN_RCPP
    Rcpp::IntegerVector totalval(_totalval);
    Rcpp::IntegerVector totallen(_totallen);
    Rcpp::NumericVector P(_P);
    Rcpp::NumericVector ambient(_ambient);

    if (totalval.size()!=totallen.size()) {
        throw std::runtime_error("length of run value/length vectors are not the same");
    }
    const int nvalues=std::accumulate(totallen.begin(), totallen.end(), 0);
    if (nvalues!=P.size()) { 
        throw std::runtime_error("sum of run lengths does not equal length of 'P' vector");
    }
   
    const int niter=check_integer_scalar(_iter, "number of iterations");
    if (niter < 1) {
        throw std::runtime_error("number of iterations should be a positive integer");
    }

    // Setting up temporary objects.
    Rcpp::IntegerVector above(nvalues);
    const size_t ngenes=ambient.size();
    if (ngenes==0) { 
        return above;
    }
    std::vector<int> tracker(ngenes);
    std::vector<double> logprob(ambient.begin(), ambient.end());
    for (auto& l : logprob) {
        l=std::log(l);
    }
    std::vector<double> cumprob(ngenes);
    std::partial_sum(ambient.begin(), ambient.end(), cumprob.begin());
    const double sumprob=cumprob.back();
    Rcpp::RNGScope rng; // after the IntegerVector!

    // Looping across iterations. 
    for (int it=0; it<niter; ++it) {
        int curtotal=0;
        std::fill(tracker.begin(), tracker.end(), 0);
        double curp=0;

        auto abIt=above.begin();
        auto tvIt=totalval.begin();
        auto tlIt=totallen.begin();
        auto pIt=P.begin();

        while (tvIt!=totalval.end()) {
            const auto& curlen=*tlIt;
            const auto& curval=*tvIt;
            
            // Sampling more points to reach the current total count.
            while (curtotal<curval) {
                auto chosen=std::lower_bound(cumprob.begin(), cumprob.end(), R::unif_rand() * sumprob) - cumprob.begin();
                if (chosen >= ngenes) { chosen=ngenes-1; }
                curp += logprob[chosen] - std::log(++tracker[chosen]);
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
    END_RCPP
}

