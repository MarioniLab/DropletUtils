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
        double& current=(*sdIt);
        ++sdIt;

        for (const auto& p : prop) { 
            double curmean=t*p;
            int out=R::rpois(curmean);
            if (out) { 
                current+=out * std::log(out/curmean) - out;
            }
            current+=curmean;
        }
    }
    
    return sum_dev;
    END_RCPP
}

SEXP calculate_pval (SEXP _obstotal, SEXP _obsresid, SEXP _simtotal, SEXP _simresid, SEXP _tol) {
    BEGIN_RCPP
    Rcpp::NumericVector obstotal(_obstotal);
    Rcpp::NumericVector obsresid(_obsresid);
    Rcpp::NumericVector simtotal(_simtotal);
    Rcpp::NumericVector simresid(_simresid);

    double last_obs=0;
    const size_t nobs=obstotal.size();
    if (nobs) { last_obs=obstotal[0]; }
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
    Rcpp::IntegerVector output(nobs);
    auto oIt=output.begin();
    Rcpp::IntegerVector totals(nobs);
    auto tIt=totals.begin();
    auto orIt=obsresid.begin();

    for (const auto& ot : obstotal) {
        const double lower=ot*tol, upper=ot/tol;
        const double& curresid=*orIt;
        ++orIt;

        // Finding the upper and lower bound.
        auto left=std::lower_bound(simtotal.begin(), simtotal.end(), lower);
        auto right=std::upper_bound(simtotal.begin(), simtotal.end(), upper);
        *tIt=right-left;
        ++tIt;

        int& nabove=(*oIt=0);
        ++oIt;
        auto srIt=simresid.begin() + (left-simtotal.begin());
        while (left!=right) {
            if (*srIt >= curresid) { ++nabove; }
            ++left;
            ++srIt;
        }
    }
         
    return Rcpp::List::create(output, totals); 
    END_RCPP
}

}
