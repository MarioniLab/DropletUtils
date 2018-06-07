#include "DropletUtils.h"

/* Defining some general-purpose downsampling functions. */

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
void downsample_counts (IN freqIt, IN freqEnd, OUT freqOut, 
        int num_total, int num_sample, 
        int& num_processed, int& num_selected) {        

    if (freqIt==freqEnd) { 
        return;
    }
    int end_of_run=num_processed + *freqIt;
    ++freqIt;

    while (num_selected < num_sample) {

        // Finding the next event with a non-zero frequency.
        while (end_of_run==num_processed && freqIt!=freqEnd) {
            end_of_run+=*freqIt;
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
void downsample_counts (IN freqIt, IN freqEnd, OUT oIt, double prop) { 
    const int num_total=std::accumulate(freqIt, freqEnd, 0), num_sample=std::round(prop*num_total);
    int num_processed=0, num_selected=0;
    downsample_counts(freqIt, freqEnd, oIt, num_total, num_sample, num_processed, num_selected);
    return;
}

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

template <class M, class O>
void downsample_matrix_internal(M mat, O output, Rcpp::NumericVector prop, Rcpp::LogicalVector bycol) {
    // Checking inputs:
    const size_t ngenes=mat->get_nrow();
    Rcpp::IntegerVector incoming(ngenes), outgoing(ngenes);
    const size_t ncells=mat->get_ncol();

    int num_total=0, num_sample=0, num_processed=0, num_selected=0;
    const bool percol=check_downsampling_mode(ncells, prop, bycol);
    if (!percol) {
        // Getting the total sum of counts in the matrix.
        for (size_t i=0; i<ncells; ++i) {
            mat->get_col(i, incoming.begin());
            num_total+=std::accumulate(incoming.begin(), incoming.end(), 0);
        }
        num_sample=std::round(num_total*prop[0]);
    }

    // Iterating across cells and downsampling the count matrix.
    auto pIt=prop.begin();
    Rcpp::RNGScope _rng; 

    for (size_t i=0; i<ncells; ++i) {
        mat->get_col(i, incoming.begin());

        // Setting up the output vector.
        std::fill(outgoing.begin(), outgoing.end(), 0);
        if (percol) { 
            downsample_counts(incoming.begin(), incoming.end(), outgoing.begin(), *pIt);
            ++pIt;
        } else {
            downsample_counts(incoming.begin(), incoming.end(), outgoing.begin(), num_total, num_sample, num_processed, num_selected);
        }

        output->set_col(i, outgoing.begin());
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
    auto otype=beachmat::output_param(rmat, false, true);
    if (rtype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(rmat);
        auto out=beachmat::create_integer_output(mat->get_nrow(), mat->get_ncol(), otype);
        downsample_matrix_internal(mat.get(), out.get(), prop, bycol);
        return out->yield();
    } else {
        auto mat=beachmat::create_numeric_matrix(rmat);
        auto out=beachmat::create_numeric_output(mat->get_nrow(), mat->get_ncol(), otype);
        downsample_matrix_internal(mat.get(), out.get(), prop, bycol);
        return out->yield();
    }
    END_RCPP    
}

/*************************************************
 ***** Downsampling (each run of) a vector. ******
 *************************************************/

SEXP downsample_runs(SEXP _groups, SEXP _reads, SEXP _prop, SEXP _bycol) {
    BEGIN_RCPP

    // Checking all of the inputs.
    Rcpp::IntegerVector groups(_groups);
    Rcpp::IntegerVector reads(_reads);
    const int nmolecules=std::accumulate(groups.begin(), groups.end(), 0);
    if (nmolecules!=reads.size()) {
        throw std::runtime_error("length of 'reads' vector should be equal to sum of RLE lengths");
    }

    int num_total=0, num_sample=0, num_processed=0, num_selected=0;
    Rcpp::NumericVector prop(_prop);
    const bool percol=check_downsampling_mode(groups.size(), prop, _bycol);
    if (!percol) {
        // Getting the total sum of counts in the vector.
        num_total=std::accumulate(reads.begin(), reads.end(), 0);
        num_sample=std::round(num_total*prop[0]);
    }

    // Setting up the output.
    Rcpp::IntegerVector output(nmolecules, 0);
    auto oIt=output.begin();
    auto rIt=reads.begin();
    auto pIt=prop.begin();

    // Iterating across the molecule groups and downsampling.
    for (const auto& g : groups) {
        if (percol) { 
            downsample_counts(rIt, rIt+g, oIt, *pIt);
            ++pIt;
        } else {
            downsample_counts(rIt, rIt+g, oIt, num_total, num_sample, num_processed, num_selected);
        }
        rIt+=g;
        oIt+=g;        
    }

    return output;
    END_RCPP;
}

