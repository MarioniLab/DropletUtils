#include "DropletUtils.h"

/* Defining some general-purpose downsampling functions. */

template<class IN, class OUT> 
void downsample_counts (IN iIt, IN iend, OUT oIt, 
        int num_total, 
        int num_sample, 
        int& offset,      // the location of 'iIt' in the total set of points.
        int& num_selected // number of points already selected.
        ) {        

    if (iIt==iend) { 
        return;
    }
    int cumulative=offset + *iIt;
    ++iIt;

    // Sampling scheme adapted from John D. Cook, https://stackoverflow.com/a/311716/15485.
    while (num_selected < num_sample) {

        /* Advancing to that point's "index" (if we had instantiated the full [0, num_total) array).
         * Note that we need to use a while loop just in case there's a whole bunch of zeroes.
         */
        while (cumulative==offset && iIt!=iend) {
            cumulative+=(*iIt);
            ++iIt;
            ++oIt;
        }

        // Breaking if all points have been iterated over.
        if (cumulative==offset && iIt==iend) { 
            break;
        }

        // Deciding whether or not to keep this point.
        if ( (num_total - offset)*R::unif_rand() < num_sample - num_selected) {
            ++(*oIt);
            ++num_selected;
        }
      
        ++offset; 
    }
    return;
}  

template<class IN, class OUT> 
void downsample_counts (IN iIt, IN iend, OUT oIt, double prop) { 
    // Convenience wrapper, when we're just downsampling in a vector.
    const int num_total=std::accumulate(iIt, iend, 0), num_sample=std::round(prop*num_total);
    int offset=0, num_selected=0;
    downsample_counts(iIt, iend, oIt, num_total, num_sample, offset, num_selected);
    return;
}

bool check_downsampling_mode (size_t ncells, Rcpp::NumericVector prop, Rcpp::LogicalVector bycol) {
    // Choosing between global downsampling or cell-specific downsampling.
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

    int num_total=0, num_sample=0, offset=0, num_selected=0;
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
            downsample_counts(incoming.begin(), incoming.end(), outgoing.begin(), num_total, num_sample, offset, num_selected);
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

    int num_total=0, num_sample=0, offset=0, num_selected=0;
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
            downsample_counts(rIt, rIt+g, oIt, num_total, num_sample, offset, num_selected);
        }
        rIt+=g;
        oIt+=g;        
    }

    return output;
    END_RCPP;
}

