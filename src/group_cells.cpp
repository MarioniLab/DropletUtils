#include "DropletUtils.h"

#include <stdexcept>
#include <algorithm>
#include <deque>

SEXP group_cells (SEXP cells, SEXP gems) {
    BEGIN_RCPP
    Rcpp::StringVector Cells(cells);
    Rcpp::IntegerVector Gems(gems);
    const size_t N=Cells.size();
    if (N!=Gems.size()) {
        throw std::runtime_error("cell and gem ID vectors should be the same length");
    }

    /* Figuring out the order (using a stable_sort for an easier comparison to downsampleMatrix in test-downsample.R).
     * Do NOT create a vector of struct's with a Rcpp::String member to improve cache locality during sorting;
     * these are really expensive to destroy, for some reason (see discussion at tidyverse/tidyselect#56).
     */
    std::vector<size_t> ordering(N);
    std::iota(ordering.begin(), ordering.end(), 0);

    std::stable_sort(ordering.begin(), ordering.end(), [&](const int& left, const int& right) {
        if (Cells[left] < Cells[right]) {
            return true;
        } else if (Cells[left] > Cells[right]) {
           return false;
        }
        return Gems[left] < Gems[right];
    });

    // Now figuring out which ones are unique.
    std::deque<size_t> unique_idx;
    std::deque<int> unique_num;

    if (N) {
        // Avoid creating lots of Rcpp::String objects, the cache freaks out and overflows.
        auto orderIt=ordering.begin();
        auto cellIt=Cells.begin()+*orderIt;
        auto gemIt=Gems.begin()+*orderIt;

        unique_idx.push_back(*orderIt);
        unique_num.push_back(1);
        ++orderIt;

        for (size_t i=1; i<N; ++i, ++orderIt) {
            auto altCellIt=Cells.begin() + *orderIt;
            auto altGemIt=Gems.begin() + *orderIt;

            if (*altCellIt!=*cellIt || *altGemIt!=*gemIt) {
                cellIt=altCellIt;
                gemIt=altGemIt;
                unique_idx.push_back(*orderIt);
                unique_num.push_back(1);
            } else {
                ++(unique_num.back());
            }
        }

        // Getting back to 1-based indexing.
        for (auto& o : ordering) { ++o; }
        for (auto& i : unique_idx) { ++i; }
    }

    Rcpp::RObject order_output, unique_output;
    if (N <= 2147483647) { // If the indices are greater than .Machine$integer.max, we switch to double.
        order_output=Rcpp::IntegerVector(ordering.begin(), ordering.end());
        unique_output=Rcpp::IntegerVector(unique_idx.begin(), unique_idx.end());
    } else {
        order_output=Rcpp::NumericVector(ordering.begin(), ordering.end());
        unique_output=Rcpp::NumericVector(unique_idx.begin(), unique_idx.end());
    }

    return Rcpp::List::create(order_output, unique_output, 
        Rcpp::IntegerVector(unique_num.begin(), unique_num.end()));
    END_RCPP
}
