#include "Rcpp.h"
#include <cmath>
#include <algorithm>

extern "C" {

SEXP calculate_random_dev (SEXP _totals, SEXP _prop) {
    BEGIN_RCPP
    Rcpp::NumericVector totals(_totals);
    Rcpp::NumericVector prop(_prop);

    const size_t ngenes=prop.size();
    const size_t npts=totals.size();
    Rcpp::NumericVector sum_dev(npts);
    Rcpp::RNGScope rng;

    auto sdIt=sum_dev.begin();
    for (const auto& t : totals) {
        double& current=(*sdIt=t); // Taking advantage of the formulation of the deviance.
        ++sdIt;

        for (const auto& p : prop) { 
            double curmean=t*p;
            int out=R::rpois(curmean);
            if (out) { current+=out * std::log(out/curmean) - out; }
        }
    }
    
    return sum_dev;
    END_RCPP
}

SEXP calculate_pval (SEXP _obstotal, SEXP _obsresid, SEXP _simtotal, SEXP _simresid, SEXP _tol) {
    BEGIN_RCPP
    Rcpp::IntegerVector obstotal(_obstotal);
    Rcpp::NumericVector obsresid(_obsresid);
    Rcpp::NumericVector simtotal(_simtotal);
    Rcpp::NumericVector simresid(_simresid);

    const size_t nobs=obstotal.size();
    if (nobs!=obsresid.size()) { throw std::runtime_error("observed vector lengths are not the same"); }
    if (simtotal.size()!=simresid.size()) { 
        throw std::runtime_error("simulated vector lengths are not the same"); 
    }
   
    Rcpp::NumericVector TOL(_tol);
    if (TOL.size()!=1) { 
        throw std::runtime_error("tolerance should be a numeric scalar");
    }
    const double tol=TOL[0];
    if (tol >= 1) {
        throw std::runtime_error("tolerance should be less than 1");
    }

    // Setting up output constructs.
    Rcpp::IntegerVector above(nobs);
    auto aIt=above.begin();
    Rcpp::IntegerVector totals(nobs);
    auto tIt=totals.begin();

    auto orIt=obsresid.begin();
    auto otIt=obstotal.begin();
    auto stLeft=simtotal.begin(), stRight=stLeft;
    auto srLeft=simresid.begin(), srRight=srLeft;

    while (otIt!=obstotal.end()) { 
        // Figuring out how many other points have similar totals.
        auto endtIt=otIt;
        while (endtIt!=obstotal.end() && *endtIt==*otIt) { 
            ++endtIt;
        }
        const size_t numO=(endtIt-otIt);
        auto endrIt=orIt + numO;

        // Finding the upper and lower bound.
        const double ot=*otIt;
        const double lower=ot*tol, upper=ot/tol;
        while (*stLeft < lower && stLeft!=simtotal.end()) { 
            ++stLeft;
            ++srLeft;
        }
        while (*stRight < upper && stRight!=simtotal.end()) { 
            ++stRight;
            ++srRight;
        }

        // Counting the number of residuals larger than the current.
        // Done with a binary search on the sorted observed residuals with the same total.
        for (auto srIt=srLeft; srIt!=srRight; ++srIt) {
            const size_t hit=std::upper_bound(orIt, endrIt, *srIt) - orIt;
            if (hit) { *(aIt+hit-1)+=1; }
        }
        if (numO>=2) { 
            // Cumulative sum to count the number of larger residuals.
            for (int i=int(numO)-2; i>=0; --i) {
                *(aIt+i)+=*(aIt+i+1);
            }
        }
        const int numS=stRight-stLeft;
        std::fill(tIt, tIt+numO, numS);
    
        // Jumping ahead.
        aIt+=numO; 
        tIt+=numO;
        orIt=endrIt;
        otIt=endtIt;
    }
         
    return Rcpp::List::create(above, totals); 
    END_RCPP
}

}
