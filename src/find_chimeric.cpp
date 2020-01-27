#include "Rcpp.h"

#include <stdexcept>
#include <vector>
#include <algorithm>

struct molecule {
    molecule (size_t i, int u) : index(i), umi(u) {}

    // need to handle situations where one sample has >2e9 UMIs.
    // otherwise, using ints for memory efficiency.
    size_t index; 
    int umi; 
};

//[[Rcpp::export(rng=false)]]
Rcpp::List find_chimeric(Rcpp::StringVector cells, Rcpp::IntegerVector umis, 
    Rcpp::IntegerVector reads, double minfrac, bool diagnostics)
{
    auto nmolecules=cells.size();
    if (nmolecules!=umis.size() || nmolecules!=reads.size()) {
        throw std::runtime_error("'reads', 'umis', 'cells' should be of the same length");
    }

    std::vector<molecule> ordering;
    ordering.reserve(nmolecules);
    auto uIt=umis.begin();
    for (size_t i=0; i<nmolecules; ++i, ++uIt) {
        ordering.push_back(molecule(i, *uIt));
    }

    // Sorting the indices.
    std::sort(ordering.begin(), ordering.end(), [&](const molecule& left, const molecule& right) {
        if (left.umi < right.umi) {
            return true;
        } else if (left.umi > right.umi) {
            return false;
        }

        // Referencing the StringVectors to avoid copying them.
        return cells[left.index] < cells[right.index];
    });

    // Setting up the output, indicating which values to keep from each sample.
    Rcpp::LogicalVector notchimeric(nmolecules);

    auto same_combination = [&] (const molecule& left, const molecule& right) {
        return left.umi==right.umi && cells[left.index]==cells[right.index];
    };

    // Iterating across runs of the same UMI/gene/cell combination.
    auto ostart=ordering.begin(), oend=ordering.begin();
    size_t nunique=0;

    while (ostart!=ordering.end()) {
        int max_nread=reads[ostart->index];
        int total_nreads=max_nread;
        auto best_mol=ostart;

        // ostart is always equal to oend at this point, so incrementing the latter to get to the next read.
        ++oend; 
        while (oend!=ordering.end() && same_combination(*ostart, *oend)) { 
            const int current_nread=reads[oend->index];
            if (current_nread > max_nread) {
                max_nread=current_nread;
                best_mol=oend;
            }

            total_nreads += current_nread;
            ++oend;
        }

        if (double(max_nread)/total_nreads >= minfrac) {
            notchimeric[best_mol->index]=1;
        }
        if (diagnostics) {
            ++nunique;
        }
        ostart=oend;
    }

    // Creating the output list.
    Rcpp::List output(2);
    output[0]=notchimeric;
    output[1]=R_NilValue;

    // Storing diagnostic information about each unique combination.
    if (diagnostics) {
        Rcpp::IntegerVector counts(nunique), totals(nunique);
        Rcpp::NumericVector prop(nunique);

        auto ostart=ordering.begin(), oend=ordering.begin();
        size_t counter=0;
        while (ostart!=ordering.end()) {
            auto& curcount=counts[counter];
            auto& curtotal=totals[counter];
            auto& curprop=prop[counter];

            while (oend!=ordering.end() && (ostart==oend || same_combination(*ostart, *oend))) { 
                const int curread=reads[oend->index];
                curtotal+=curread;
                if (curread > curprop) { curprop=curread; }
                ++curcount;
                ++oend;
            }

            if (curtotal) {
                curprop/=curtotal;
            }

            ostart=oend;
            ++counter;
        }
        output[1]=Rcpp::List::create(counts, totals, prop);
    }

    return output;
}
