#include "Rcpp.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/utils/const_column.h"
#include "utils.h"

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <cstdint>

// Using big integers to avoid overflow when summing across a large matrix.

typedef std::uint64_t bigint_t;

template<class IT>
bigint_t bigsum(IT start, IT end) {
    bigint_t sum=0;
    while (start!=end) {
        sum+=static_cast<int>(*start);
        ++start;
    }
    return sum;
}

/* Defining some general-purpose downsampling functions. */

class downsampler {
public:
    void set_global(bigint_t total, double prop) {
        num_total=total;
        set_num_sample(prop);
    }

    /* This class considers sampling events without replacement from a vector.
     * Here, though, the vector contains frequencies of events rather than the events themselves.
     * The sampling scheme is adapted from John D. Cook, https://stackoverflow.com/a/311716/15485.
     * 
     * freqIt: An iterator pointing to the start of the frequency vector.
     * freqEnd: An interator pointing to the end of the frequency vector.
     * freqOut: An iterator pointing to an output vector, indicating how many instances of each event have been sampled.
     * num_total: An integer scalar specifying the total number of all events.
     * num_sample: An integer scalar specifying the number of events to sample without replacement.
     * num_processed: An integer scalar specifying the number of events that have already been considered for selection.
     * num_selected: An integer scalar specifying the number of events that have already been selected.
     * 
     * Note that num_total may not be simply a sum of all values from [freqIt, freqEnd).
     * This is because we allow multiple applications of this function to sample without replacement from a series of vectors.
     * We keep track of 'num_processed' and 'num_selected' to ensure correct sampling when moving from one vector to another.
     */
    template<class IN, class OUT> 
    void operator()(IN freqIt, IN freqEnd, OUT freqOut) {        
        while (freqIt!=freqEnd && num_selected < num_sample) {
            for (int i=0; i<*freqIt && num_sample > num_selected; ++i) {
                // Deciding whether or not to keep this instance of this event.
                // This is a safe way of computing NUM_YET_TO_SELECT/NUM_YET_TO_PROCESS > runif(1), 
                // avoiding issues with integer division.
                if ( (num_total - num_processed)*R::unif_rand() < num_sample - num_selected) {
                    ++(*freqOut);
                    ++num_selected;
                }
                ++num_processed;
            }
         
            ++freqIt;
            ++freqOut;
        }
        return;
    }  

    /* Convenience wrapper when we're just downsampling in a single vector.
     * In this case, num_total is just a sum of [freqIt, freqEnd). 
     * There is also no need for any special values of num_processed and 'num_selected'.
     */
    template<class IN, class OUT> 
    void operator()(IN freqIt, IN freqEnd, OUT oIt, double prop) { 
        num_total=bigsum(freqIt, freqEnd);
        set_num_sample(prop);
        num_processed=0;
        num_selected=0;

        (*this)(freqIt, freqEnd, oIt);
        return;
    }

private:
    bigint_t num_total=0, num_sample=0, num_processed=0, num_selected=0;

    void set_num_sample(double prop) {
        num_sample=std::round(prop*num_total);
    }
};

void check_downsampling_mode (size_t ncells, Rcpp::NumericVector prop, bool bycol) 
// Choosing between global downsampling or cell-specific downsampling.
{
    if (bycol) { 
        if (prop.size()!=ncells) {
            throw std::runtime_error("length of 'prop' should be equal to number of cells when 'bycol=TRUE'");
        }
        for (const auto& curprop : prop) { 
            if (curprop < 0 || curprop > 1) { 
                throw std::runtime_error("downsampling proportion must lie in [0, 1]");
            }
        }
    } else {
        if (prop.size()!=1) {
            throw std::runtime_error("downsampling proportion should be a numeric scalar when 'bycol=FALSE'");
        }
        double curprop=prop[0];
        if (curprop < 0 || curprop > 1) { 
            throw std::runtime_error("downsampling proportion must lie in [0, 1]");
        }
    }
}

/*************************************************
 **** Downsampling (each column of) a matrix. ****
 *************************************************/

template <typename V, class M, class O>
Rcpp::RObject downsample_matrix_internal(Rcpp::RObject input, Rcpp::NumericVector prop, bool bycol) {
    auto mat=beachmat::create_matrix<M>(input);
    auto otype=beachmat::output_param(mat.get());
    auto output=beachmat::create_output<O>(mat->get_nrow(), mat->get_ncol(), otype);

    // Checking inputs:
    const size_t ngenes=mat->get_nrow();
    const size_t ncells=mat->get_ncol();
    V incoming(ngenes);
    Rcpp::IntegerVector outgoing(ngenes);

    beachmat::const_column<M> col_holder(mat.get());

    // Setting variables for global downsampling, if desired.
    check_downsampling_mode(ncells, prop, bycol);

    downsampler down;
    if (!bycol) {
        bigint_t num_total=0;
        for (size_t i=0; i<ncells; ++i) {
            col_holder.fill(i);
            auto it=col_holder.get_values();
            num_total+=bigsum(it, it+col_holder.get_n());
        }
        down.set_global(num_total, prop[0]);
    }

    /* Iterating across cells and downsampling the count matrix.
     * Note that the RNGscope destructor may trigger a garbage collection,
     * so we enclose it in its own scope to ensure that it doesn't 
     * collect the unprotected output of yield().
     */
    {
        auto pIt=prop.begin();
        Rcpp::RNGScope _rng; 

        for (size_t i=0; i<ncells; ++i) {
            col_holder.fill(i);
            auto valS = col_holder.get_values();
            auto valE = valS + col_holder.get_n();

            // Downsampling.
            if (bycol) { 
                down(valS, valE, outgoing.begin(), *pIt);
                ++pIt;
            } else {
                down(valS, valE, outgoing.begin());
            }

            // Saving and then clearing the output vector.
            if (col_holder.is_sparse()) {
                output->set_col_indexed(i, col_holder.get_n(), col_holder.get_indices(), outgoing.begin());
                std::fill(outgoing.begin(), outgoing.begin() + col_holder.get_n(), 0);
            } else {
                output->set_col(i, outgoing.begin());
                std::fill(outgoing.begin(), outgoing.end(), 0);
            }
        }
    }

    return output->yield();
}

//[[Rcpp::export]]
Rcpp::RObject downsample_matrix(Rcpp::RObject rmat, Rcpp::NumericVector prop, bool bycol) {
    int rtype=beachmat::find_sexp_type(rmat);
    if (rtype==INTSXP) {
        return downsample_matrix_internal<Rcpp::IntegerVector, 
           beachmat::integer_matrix, 
           beachmat::integer_output>(rmat, prop, bycol);
    } else {
        return downsample_matrix_internal<Rcpp::NumericVector, 
           beachmat::numeric_matrix, 
           beachmat::numeric_output>(rmat, prop, bycol);
    }
}

/*************************************************
 ***** Downsampling (each run of) a vector. ******
 *************************************************/

//[[Rcpp::export]]
Rcpp::IntegerVector downsample_runs(Rcpp::IntegerVector cells, Rcpp::IntegerVector reads, 
    Rcpp::NumericVector prop, bool bycol) 
{
    // Checking all of the inputs.
    const bigint_t nmolecules=bigsum(cells.begin(), cells.end());
    if (nmolecules!=reads.size()) {
        throw std::runtime_error("length of 'reads' vector should be equal to sum of RLE lengths");
    }

    downsampler down;
    check_downsampling_mode(cells.size(), prop, bycol);
    if (!bycol) {
        down.set_global(bigsum(reads.begin(), reads.end()), prop[0]);
    }

    // Setting up the output.
    Rcpp::IntegerVector output(nmolecules);
    auto oIt=output.begin();
    auto rIt=reads.begin();
    auto pIt=prop.begin();

    // Iterating across the molecule cells and downsampling.
    for (const auto& cell : cells) {
        if (bycol) { 
            down(rIt, rIt+cell, oIt, *pIt);
            ++pIt;
        } else {
            down(rIt, rIt+cell, oIt);
        }
        rIt+=cell;
        oIt+=cell;
    }

    return output;
}
