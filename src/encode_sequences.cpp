#include "Rcpp.h"
#include <stdexcept>
#include <sstream>

//[[Rcpp::export(rng=false)]]
Rcpp::IntegerVector encode_sequences (Rcpp::StringVector Seqs) {
    Rcpp::IntegerVector output(Seqs.size());
    for (size_t i=0; i<output.size(); ++i) {
        Rcpp::String current=Seqs[i];
        int& curenc=output[i];

        const char * ptr=current.get_cstring();
        const size_t len=Rf_length(current.get_sexp());
        if (len > 15) {
            throw std::runtime_error("32-bit integers do not support sequences longer than 15 nt");
        }

        size_t mult=1;
        for (size_t j=0; j<len; ++j) {
            switch (ptr[len-j-1]) {
                case 'A':
                    break;
                case 'C':
                    curenc+=mult;
                    break;
                case 'G':
                    curenc+=mult*2;
                    break;
                case 'T':
                    curenc+=mult*3;
                    break;
                default:
                    {
                        std::stringstream err;
                        err << "unrecognized character in '" << ptr << "'";
                        throw std::runtime_error(err.str());
                    }
            }
            mult*=4;
        }
    }

    return output;
}
