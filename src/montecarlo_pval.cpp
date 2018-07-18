#include "DropletUtils.h"

SEXP montecarlo_pval (SEXP totalval, SEXP totallen, SEXP prob, SEXP ambient, SEXP iter) {
    BEGIN_RCPP
    Rcpp::IntegerVector Totalval(totalval);
    Rcpp::IntegerVector Totallen(totallen);
    Rcpp::NumericVector Prob(prob);
    Rcpp::NumericVector Ambient(ambient);

    if (Totalval.size()!=Totallen.size()) {
        throw std::runtime_error("length of run value/length vectors are not the same");
    }
    const int nvalues=std::accumulate(Totallen.begin(), Totallen.end(), 0);
    if (nvalues!=Prob.size()) {
        throw std::runtime_error("sum of run lengths does not equal length of 'P' vector");
    }

    const int niter=check_integer_scalar(iter, "number of iterations");
    if (niter < 1) {
        throw std::runtime_error("number of iterations should be a positive integer");
    }

    // Setting up temporary objects.
    Rcpp::IntegerVector above(nvalues);
    const size_t ngenes=Ambient.size();
    if (ngenes==0) {
        return above;
    }

    std::vector<int> tracker(ngenes);
    std::vector<double> logprob(Ambient.begin(), Ambient.end());
    for (auto& l : logprob) {
        l=std::log(l);
    }
    std::vector<double> cumprob(ngenes);
    std::partial_sum(Ambient.begin(), Ambient.end(), cumprob.begin());
    const double sumprob=cumprob.back();
    Rcpp::RNGScope rng; // after the IntegerVector!

    // Looping across iterations.
    for (int it=0; it<niter; ++it) {
        int curtotal=0;
        std::fill(tracker.begin(), tracker.end(), 0);
        double curp=0;

        auto abIt=above.begin();
        auto tvIt=Totalval.begin();
        auto tlIt=Totallen.begin();
        auto pIt=Prob.begin();

        while (tvIt!=Totalval.end()) {
            const auto& curlen=*tlIt;
            const auto& curval=*tvIt;

            // Sampling more points to reach the current total count.
            while (curtotal<curval) {
                auto chosen=std::lower_bound(cumprob.begin(), cumprob.end(), R::unif_rand() * sumprob) - cumprob.begin();
                if (chosen >= ngenes) {
                    chosen=ngenes-1; // Some protection against the very-unlikely case.
                }
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
    for (const auto& curlen : Totallen) {
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

SEXP montecarlo_pval_alpha (SEXP totalval, SEXP totallen, SEXP prob, SEXP ambient, SEXP iter, SEXP alpha) {
    BEGIN_RCPP
    Rcpp::IntegerVector Totalval(totalval);
    Rcpp::IntegerVector Totallen(totallen);
    Rcpp::NumericVector Prob(prob);
    Rcpp::NumericVector Ambient(ambient);

    if (Totalval.size()!=Totallen.size()) {
        throw std::runtime_error("length of run value/length vectors are not the same");
    }
    const int nvalues=std::accumulate(Totallen.begin(), Totallen.end(), 0);
    if (nvalues!=Prob.size()) {
        throw std::runtime_error("sum of run lengths does not equal length of 'P' vector");
    }

    const int niter=check_integer_scalar(iter, "number of iterations");
    if (niter < 1) {
        throw std::runtime_error("number of iterations should be a positive integer");
    }

    double Alpha=check_numeric_scalar(alpha, "alpha");
    const bool use_alpha=R_FINITE(Alpha);
    if (use_alpha && Alpha <= 0) {
        throw std::runtime_error("alpha must be positive if specified");
    }

    // Setting up temporary objects.
    Rcpp::IntegerVector above(nvalues);
    const size_t ngenes=Ambient.size();
    if (ngenes==0) {
        return above;
    }

    std::vector<int> tracker(ngenes);
    std::vector<double> tmpprob(ngenes), cumprob(ngenes);
    Rcpp::RNGScope rng; // after the IntegerVector!

    // Looping across iterations, using a new probability vector per iteration.
    for (int it=0; it<niter; ++it) {
        for (size_t ldx=0; ldx<ngenes; ++ldx) {
            tmpprob[ldx]=R::rgamma(Ambient[ldx] * Alpha, 1);
        }
        std::partial_sum(tmpprob.begin(), tmpprob.end(), cumprob.begin());
        const double sumprob=cumprob.back();

        int curtotal=0;
        std::fill(tracker.begin(), tracker.end(), 0);
        double curp=0;

        auto abIt=above.begin();
        auto tvIt=Totalval.begin();
        auto tlIt=Totallen.begin();
        auto pIt=Prob.begin();

        while (tvIt!=Totalval.end()) {
            const auto& curlen=*tlIt;
            const auto& curval=*tvIt;

            // Sampling more points to reach the current total count.
            while (curtotal<curval) {
                auto chosen=std::lower_bound(cumprob.begin(), cumprob.end(), R::unif_rand() * sumprob) - cumprob.begin();
                if (chosen >= ngenes) {
                    chosen=ngenes-1; // Some protection against the very-unlikely case.
                }

                auto& curnum=tracker[chosen];
                curp += std::log(Ambient[chosen] * Alpha + curnum); // Difference of Gamma's in data-dependent component of Dirichlet-multinomial.
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
    for (const auto& curlen : Totallen) {
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

