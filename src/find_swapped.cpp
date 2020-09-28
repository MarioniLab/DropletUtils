#include "Rcpp.h"
#include "beachmat3/beachmat.h"
#include "utils.h"

#include <stdexcept>
#include <vector>
#include <deque>
#include <algorithm>

template<class V>
std::vector<V> process_list(Rcpp::List incoming) {
    const size_t nsamples=incoming.size();
    std::vector<V> output(nsamples);
    for (size_t i=0; i<output.size(); ++i) {
        output[i] = V(incoming[i]);
    }
    return(output);
}

template<class U, class V>
void compare_lists(U left, V right) {
    if (left.size()!=right.size()) {
        throw std::runtime_error("lists are not of the same length");
    }
    const size_t nsamples=left.size();
    for (size_t i=0; i<nsamples; ++i) {
        if (left[i].size()!=right[i].size()) {
            throw std::runtime_error("list vectors are not of the same length");
        }
    }
    return;
}

struct molecule {
    molecule (int s, size_t i, int g, int u) : index(i), sample(s), gene(g), umi(u) {}

    // need to handle situations where one sample has >2e9 UMIs.
    // otherwise, using ints for memory efficiency.
    size_t index; 
    int sample, gene, umi; 
};

/* Identifies which molecules should be retained in which samples,
 * given the cell, gene and UMI combination for each molecule per sample.
 * Also returns a diagnostic matrix of molecule-sample read counts.
 */

//[[Rcpp::export(rng=false)]]
Rcpp::List find_swapped(Rcpp::List cells, Rcpp::List genes, Rcpp::List umis, Rcpp::List reads, double minfrac, bool diagnostics) {
    auto Cells=process_list<Rcpp::StringVector>(cells);
    auto Genes=process_list<Rcpp::IntegerVector>(genes);
    auto Umis=process_list<Rcpp::IntegerVector>(umis);
    auto Reads=process_list<Rcpp::IntegerVector>(reads);

    compare_lists(Cells, Genes);
    compare_lists(Cells, Umis);
    compare_lists(Cells, Reads);

    // Setting up the ordering vector.
    const size_t nsamples=Cells.size();
    size_t nmolecules=0;
    for (size_t i=0; i<nsamples; ++i) {
        nmolecules+=Cells[i].size();
    } 

    std::vector<molecule> ordering;
    ordering.reserve(nmolecules);
    for (size_t i=0; i<nsamples; ++i) {
        const size_t cur_nmol=Cells[i].size();
        const auto& cur_genes=Genes[i];
        const auto& cur_umis=Umis[i];

        auto gIt=cur_genes.begin();
        auto uIt=cur_umis.begin();
        for (size_t j=0; j<cur_nmol; ++j, ++gIt, ++uIt) {
            ordering.push_back(molecule(i, j, *gIt, *uIt));
        }
    }
    
    // Sorting the indices.
    std::sort(ordering.begin(), ordering.end(), [&](const molecule& left, const molecule& right) {
        if (left.gene < right.gene) {
            return true;
        } else if (left.gene > right.gene) {
            return false;
        }

        if (left.umi < right.umi) {
            return true;
        } else if (left.umi > right.umi) {
            return false;
        }

        // Referencing the StringVectors to avoid copying them.
        return Cells[left.sample][left.index] < Cells[right.sample][right.index];
    });
    
    // Setting up the output, indicating which values to keep from each sample.
    std::vector<Rcpp::LogicalVector> notswapped(nsamples);
    for (size_t i=0; i<nsamples; ++i) {
        notswapped[i]=Rcpp::LogicalVector(Cells[i].size());
    }

    auto same_combination = [&] (const molecule& left, const molecule& right) {
        return left.gene==right.gene && left.umi==right.umi
            && Cells[left.sample][left.index]==Cells[right.sample][right.index];
    };

    // Iterating across runs of the same UMI/gene/cell combination.
    auto ostart=ordering.begin(), oend=ordering.begin();

    while (ostart!=ordering.end()) {
        int max_nread=Reads[ostart->sample][ostart->index];
        int total_nreads=max_nread;
        auto best_mol=ostart;

        // ostart is always equal to oend at this point, so incrementing the latter to get to the next read.
        ++oend; 
        while (oend!=ordering.end() && same_combination(*ostart, *oend)) { 
            const int current_nread=Reads[oend->sample][oend->index];
            if (current_nread > max_nread) {
                max_nread=current_nread;
                best_mol=oend;
            }

            total_nreads += current_nread;
            ++oend;
        }

        if (double(max_nread)/total_nreads >= minfrac) {
            notswapped[best_mol->sample][best_mol->index]=1;
        }
        ostart=oend;
    }

    // Creating the output list.
    Rcpp::List outlist(nsamples);
    for (size_t i=0; i<nsamples; ++i) {
        outlist[i]=notswapped[i];
    }
    Rcpp::List output(2);
    output[0]=outlist;
    output[1]=R_NilValue;

    // Storing diagnostic information about each unique combination.
    if (diagnostics) {
        std::deque<std::pair<std::pair<int, int>, double> > storage;

        auto ostart=ordering.begin(), oend=ordering.begin();
        size_t counter=0;
        while (ostart!=ordering.end()) {
            size_t nnzero=0;

            // Adding read counts per molecule, storing them in the matrix, then wiping them.
            while (oend!=ordering.end() && (ostart==oend || same_combination(*ostart, *oend))) { 
                if (nnzero >= nsamples) {
                    throw std::runtime_error("multiple instances of the same combination observed in a single sample");
                }
                storage.push_back(std::make_pair(std::make_pair(oend->sample, counter), Reads[oend->sample][oend->index]));
                ++oend;
            }

            ostart=oend;
            ++counter;
        }
        
        std::sort(storage.begin(), storage.end());
        output[1]=beachmat::as_gCMatrix<Rcpp::NumericVector>(counter, nsamples, storage);
    }

    return output;
}
