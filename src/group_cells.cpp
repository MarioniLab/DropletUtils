#include "DropletUtils.h"

SEXP group_cells (SEXP cells, SEXP gems) {
    BEGIN_RCPP
    Rcpp::StringVector Cells(cells);
    Rcpp::IntegerVector Gems(gems);
    if (Cells.size()!=Gems.size()) {
        throw std::runtime_error("cell and gem ID vectors should be the same length");
    }

    // Figuring out the order (using a stable_sort for an easier comparison to downsampleMatrix in test-downsample.R).
    Rcpp::IntegerVector output(Cells.size());
    std::iota(output.begin(), output.end(), 0);

    std::stable_sort(output.begin(), output.end(), [&](const int& left, const int& right) {
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

    if (output.size()) {
        unique_Cells.push_back(Cells[output[0]]);
        unique_Gems.push_back(Gems[output[0]]);
        unique_num.push_back(1);

        for (size_t i=1; i<output.size(); ++i) {
            auto cur_cell=Cells[output[i]];
            auto cur_gem=Gems[output[i]];
            if (cur_cell!=unique_Cells.back() || cur_gem!=unique_Gems.back()) {
                unique_Cells.push_back(cur_cell);
                unique_Gems.push_back(cur_gem);
                unique_num.push_back(1);
            } else {
                ++(unique_num.back());
            }
        }

        // Getting back to 1-based indexing.
        for (auto& o : output) { ++o; }
    }

    return Rcpp::List::create(output, R_NilValue,
        Rcpp::StringVector(unique_Cells.begin(), unique_Cells.end()),
        Rcpp::IntegerVector(unique_Gems.begin(), unique_Gems.end()),
        Rcpp::IntegerVector(unique_num.begin(), unique_num.end()));
    END_RCPP
}
