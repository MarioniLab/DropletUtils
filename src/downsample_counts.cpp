#include "DropletUtils.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
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

    /* This function considers sampling events without replacement from a vector.
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
        if (freqIt==freqEnd) { 
            return;
        }
        bigint_t end_of_run=num_processed + static_cast<int>(*freqIt);
        ++freqIt;

        while (num_selected < num_sample) {

            // Finding the next event with a non-zero frequency.
            while (end_of_run==num_processed && freqIt!=freqEnd) {
                end_of_run+=static_cast<int>(*freqIt);
                ++freqIt;
                ++freqOut;
            }

            // Breaking if all points have been iterated over.
            if (end_of_run==num_processed && freqIt==freqEnd) { 
                break;
            }

            // Deciding whether or not to keep this instance of this event.
            // This is a safe way of computing NUM_YET_TO_SELECT/NUM_YET_TO_PROCESS > runif(1), avoiding issues with integer division.
            if ( (num_total - num_processed)*R::unif_rand() < num_sample - num_selected) {
                ++(*freqOut);
                ++num_selected;
            }
         
            // Moving onto the next instance of the same event.
            ++num_processed; 
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

bool check_downsampling_mode (size_t ncells, Rcpp::NumericVector prop, Rcpp::LogicalVector bycol) 
// Choosing between global downsampling or cell-specific downsampling.
{
    const bool do_bycol=check_logical_scalar(bycol, "per-column specifier");
    if (do_bycol) { 
        if (prop.size()!=ncells) {
            throw std::runtime_error("length of 'prop' should be equal to number of cells");
        }
        for (const auto& curprop : prop) { 
            if (curprop < 0 || curprop > 1) { 
                throw std::runtime_error("downsampling proportion must lie in [0, 1]");
            }
        }
    } else {
        const double curprop=check_numeric_scalar(prop, "downsampling proportion");
        if (curprop < 0 || curprop > 1) { 
            throw std::runtime_error("downsampling proportion must lie in [0, 1]");
        }
    }
    return do_bycol;
}

/*************************************************
 **** Downsampling (each column of) a matrix. ****
 *************************************************/

template <typename V, class M, class O>
void downsample_matrix_internal(M mat, O output, Rcpp::NumericVector prop, Rcpp::LogicalVector bycol) {
    // Checking inputs:
    const size_t ngenes=mat->get_nrow();
    const size_t ncells=mat->get_ncol();
    V incoming(ngenes);
    Rcpp::IntegerVector outgoing(ngenes);

    auto raw_type=mat->col_raw_type();
    auto raws=mat->set_up_raw();
    const bool is_sparse=raw_type=="sparse";
    const bool is_dense=raw_type=="sparse";

    // Setting variables for global downsampling, if desired.
    downsampler down;
    const bool percol=check_downsampling_mode(ncells, prop, bycol);
    if (!percol) {
        bigint_t num_total=0;
        if (is_sparse) {
            for (size_t i=0; i<ncells; ++i) {
                mat->get_col_raw(i, raws);
                auto it=raws.get_values_start();
                num_total+=bigsum(it, it+raws.get_n());
            }
        } else if (is_dense) {
            for (size_t i=0; i<ncells; ++i) {
                mat->get_col_raw(i, raws);
                auto it=raws.get_values_start();
                num_total+=bigsum(it, it+ngenes);
            }
        } else {
            for (size_t i=0; i<ncells; ++i) {
                mat->get_col(i, incoming.begin());
                num_total+=bigsum(incoming.begin(), incoming.end());
            }
        }
        down.set_global(num_total, prop[0]);
    }

    // Iterating across cells and downsampling the count matrix.
    auto pIt=prop.begin();
    Rcpp::RNGScope _rng; 

    for (size_t i=0; i<ncells; ++i) {
        typename V::iterator valS, valE;

        if (is_sparse) {
            mat->get_col_raw(i, raws);
            valS=raws.get_values_start();
            valE=valS + raws.get_n();
        } else if (is_dense) {
            mat->get_col_raw(i, raws);
            valS=raws.get_values_start();
            valE=valS+ngenes;
        } else {
            valS=incoming.begin();
            mat->get_col(i, valS);
            valE=incoming.end();
        }

        // Downsampling.
        if (percol) { 
            down(valS, valE, outgoing.begin(), *pIt);
            ++pIt;
        } else {
            down(valS, valE, outgoing.begin());
        }

        // Saving and then clearing the output vector.
        if (is_sparse) {
            output->set_col_indexed(i, raws.get_n(), raws.get_structure_start(), outgoing.begin());
            std::fill(outgoing.begin(), outgoing.begin()+raws.get_n(), 0);
        } else {
            output->set_col(i, outgoing.begin());
            std::fill(outgoing.begin(), outgoing.end(), 0);
        }
    }

    /* Note that the RNGscope destructor may trigger a garbage collection.
     * I'm not sure that the object from output->yield() remains protected if a compiler does not implement RVO.
     * This would result in a copy and destruction, and a point in time at which the output memory is unprotected.
     * Hence, we do not perform an output->yield() to return out of this function.
     */
    return;
}

SEXP downsample_matrix(SEXP rmat, SEXP prop, SEXP bycol) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(rmat);
    if (rtype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(rmat);
        auto otype=beachmat::output_param(mat->get_class(), mat->get_package());
        auto out=beachmat::create_integer_output(mat->get_nrow(), mat->get_ncol(), otype);
        downsample_matrix_internal<Rcpp::IntegerVector>(mat.get(), out.get(), prop, bycol);
        return out->yield();
    } else {
        auto mat=beachmat::create_numeric_matrix(rmat);
        auto otype=beachmat::output_param(mat->get_class(), mat->get_package());
        auto out=beachmat::create_numeric_output(mat->get_nrow(), mat->get_ncol(), otype);
        downsample_matrix_internal<Rcpp::NumericVector>(mat.get(), out.get(), prop, bycol);
        return out->yield();
    }
    END_RCPP    
}

/*************************************************
 ***** Downsampling (each run of) a vector. ******
 *************************************************/

SEXP downsample_runs(SEXP cells, SEXP reads, SEXP prop, SEXP bycol) {
    BEGIN_RCPP

    // Checking all of the inputs.
    Rcpp::IntegerVector cell_vec(cells);
    Rcpp::IntegerVector read_vec(reads);
    const bigint_t nmolecules=bigsum(cell_vec.begin(), cell_vec.end());
    if (nmolecules!=read_vec.size()) {
        throw std::runtime_error("length of 'reads' vector should be equal to sum of RLE lengths");
    }

    downsampler down;
    Rcpp::NumericVector proportions(prop);
    const bool percol=check_downsampling_mode(cell_vec.size(), proportions, bycol);
    if (!percol) {
        down.set_global(bigsum(read_vec.begin(), read_vec.end()), proportions[0]);
    }

    // Setting up the output.
    Rcpp::IntegerVector output(nmolecules);
    auto oIt=output.begin();
    auto rIt=read_vec.begin();
    auto pIt=proportions.begin();

    // Iterating across the molecule cell_vec and downsampling.
    for (const auto& cell : cell_vec) {
        if (percol) { 
            down(rIt, rIt+cell, oIt, *pIt);
            ++pIt;
        } else {
            down(rIt, rIt+cell, oIt);
        }
        rIt+=cell;
        oIt+=cell;
    }

    return output;
    END_RCPP;
}

