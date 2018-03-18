#include "DropletUtils.h"

SEXP find_swapped(SEXP _groups, SEXP _reads, SEXP _minfrac) {
    BEGIN_RCPP

    // Checking all of the inputs.
    Rcpp::IntegerVector groups(_groups);
    Rcpp::IntegerVector reads(_reads);
    const int nmolecules=std::accumulate(groups.begin(), groups.end(), 0);
    if (nmolecules!=reads.size()) {
        throw std::runtime_error("length of 'reads' vector should be equal to sum of RLE lengths");
    }
    const double mf=check_numeric_scalar(_minfrac, "minimum fraction");

    // Setting up the output.
    Rcpp::LogicalVector output(nmolecules, 1);
    auto oIt=output.begin();

    // Iterating across the molecule groups.
    auto rIt=reads.begin();
    for (const auto& g : groups) {
        int topindex=0, topreads=*rIt, totalreads=topreads;
        ++rIt;

        for (int i=1; i<g; ++i, ++rIt) {
            if (topreads < *rIt) {
                topreads=*rIt;
                topindex=i;
            }
            totalreads+=*rIt;
        }
        
        if (double(totalreads)*mf <= double(topreads)) {
            *(oIt+topindex)=0;
        }
        oIt+=g;
    }

    return output;
    END_RCPP;
}
