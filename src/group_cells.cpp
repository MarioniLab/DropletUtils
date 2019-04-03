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
    std::deque<Rcpp::String> unique_Cells;
    std::deque<int> unique_Gems;
    std::deque<int> unique_num;

    if (N) {
        unique_Cells.push_back(Cells[ordering[0]]);
        unique_Gems.push_back(Gems[ordering[0]]);
        unique_num.push_back(1);

        for (size_t i=1; i<N; ++i) {
            auto cur_cell=Cells[ordering[i]];
            auto cur_gem=Gems[ordering[i]];
            if (cur_cell!=unique_Cells.back() || cur_gem!=unique_Gems.back()) {
                unique_Cells.push_back(cur_cell);
                unique_Gems.push_back(cur_gem);
                unique_num.push_back(1);
            } else {
                ++(unique_num.back());
            }
        }

        // Getting back to 1-based indexing.
        for (auto& o : ordering) { ++o; }
    }

    Rcpp::RObject output;
    if (N <= 2147483647) { // If the indices are greater than .Machine$integer.max, we switch to double.
        output=Rcpp::IntegerVector(ordering.begin(), ordering.end());
    } else {
        output=Rcpp::NumericVector(ordering.begin(), ordering.end());
    }

    return Rcpp::List::create(output, R_NilValue,
        Rcpp::StringVector(unique_Cells.begin(), unique_Cells.end()),
        Rcpp::IntegerVector(unique_Gems.begin(), unique_Gems.end()),
        Rcpp::IntegerVector(unique_num.begin(), unique_num.end()));
    END_RCPP
}
