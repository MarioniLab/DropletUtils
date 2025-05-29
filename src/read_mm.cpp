#include "eminem/eminem.hpp"
#include "byteme/byteme.hpp"
#include "subpar/subpar.hpp"

#include "Rcpp.h"

#include <vector>
#include <stdexcept>
#include <limits>

template<typename Type_>
void sort_SVT_SparseMatrix_columns(const std::vector<int*>& iptrs, const std::vector<Type_*>& vptrs, const std::vector<R_xlen_t>& num, int threads) {
    auto NC = iptrs.size();
    subpar::parallelize_range(threads, NC, [&](int, decltype(NC) start, decltype(NC) length) -> void {
        std::vector<std::pair<int, Type_> > sortbuffer;
        for (decltype(start) c = start, end = start + length; c < end; ++c) {
            auto iptr = iptrs[c];
            auto n = used[c];
            if (std::is_sorted(iptr, iptr + n)) {
                continue;
            }
            auto vptr = vptrs[c];
            sortbuffer.clear();
            for (decltype(n) i = 0; i < n; ++i) {
                sortbuffer.emplace_back(iptr[i], vptr[i]);
            }
            std::sort(sortbuffer.begin(), sortbuffer.end());
            for (decltype(n) i = 0; i < n; ++i) {
                const auto& current = sortbuffer[i];
                iptr[i] = current.first;
                vptr[i] = current.second;
            }
        }
    });
}

Rcpp::RObject read_mm_two_pass_CsparseMatrix(const std::string& path, const std::vector<R_xlen_t>& nnz_per_col, int threads) {
    auto NC = nnz_per_col.size();
    Rcpp::List contents(NC);
    std::vector<int*> iptrs(NC);
    std::vector<R_xlen_t> used(NC);

    byteme::ParseSomeFileOptions opt;
    opt.threads = threads;
    auto parser = parse_some_file(path.c_str(), opt);
    const auto& banner = parser.get_preamble();

    auto btype = banner.type();
    if (btype == eminem::Field::REAL || btype == eminem::Field::DOUBLE) {
        std::vector<double*> vptrs(NC);
        for (decltype(NC) c = 0; c < NC; ++c) {
            Rcpp::NumericVector values(nnz_per_col[c]);
            Rcpp::IntegerVector indices(nnz_per_col[c]);
            iptrs[c] = indices.begin(); // these pointers should still be valid after the std::move as they refer to R-managed allocations.
            vptrs[c] = values.begin();
            contents[c] = Rcpp::List::create(std::move(indices), std::move(values));
        }

        scan_real([&](eminem::Index r, eminem::Index c, double val) -> void {
            auto& pos = used[c];
            iptrs[c][pos] = r;
            vptrs[c][pos] = val;
            ++pos;
        });

        sort_SVT_SparseMatrix_columns(iptrs, vptrs, used, threads);

    } else if (btype == eminem::Field::INTEGER) {
        std::vector<int*> vptrs(NC);
        for (decltype(NC) c = 0; c < NC; ++c) {
            Rcpp::IntegerVector values(nnz_per_col[c]);
            Rcpp::IntegerVector indices(nnz_per_col[c]);
            iptrs[c] = indices.begin(); // these pointers should still be valid after the std::move as they refer to R-managed allocations.
            vptrs[c] = values.begin();
            contents[c] = Rcpp::List::create(std::move(indices), std::move(values));
        }

        scan_integer([&](eminem::Index r, eminem::Index c, double val) -> void {
            auto& pos = used[c];
            iptrs[c][pos] = r;
            vptrs[c][pos] = val;
            ++pos;
        });

        sort_SVT_SparseMatrix_columns(iptrs, vptrs, used, threads);

    } else {
        for (decltype(NC) c = 0; c < NC; ++c) {
            Rcpp::IntegerVector indices(nnz_per_col[c]);
            iptrs[c] = indices.begin(); // these pointers should still be valid after the std::move as they refer to R-managed allocations.
            contents[c] = Rcpp::List::create(std::move(indices), R_NilValue);
        }

        subpar::parallelize_range(threads, NC, [&](int, decltype(NC) start, decltype(NC) length) -> void {
            for (decltype(start) c = start, end = start + length; c < end; ++c) {
                auto iptr = iptrs[c];
                auto n = used[c];
                if (std::is_sorted(iptr, iptr + n)) {
                    continue;
                }
                std::sort(iptr, iptr + n);
            }
        });
    }

    return contents;
}

Rcpp::RObject read_mm_two_pass_CsparseMatrix(const std::string& path, const std::vector<R_xlen_t>& nnz_per_col, int threads) {
    auto NC = nnz_per_col.size();
    std::vector<R_xlen_t> offsets(NC + 1);
    for (decltype(NC) c = 0; c < NC; ++c) {
        offsets[c + 1] = offsets[c] + nnz_per_col[c];
    }

    auto ntotal = offsets.back();
    if (ntotal > std::numeric_limits<int>::max()) {
        throw std::runtime_error("more non-zero elements than is supported in a CsparseMatrix");
    }

    Rcpp::IntegerVector indptr(offsets.begin(), offsets.end());
    Rcpp::IntegerVector row_indices(ntotal);

    byteme::ParseSomeFileOptions opt;
    opt.threads = threads;
    auto parser = parse_some_file(path.c_str(), opt);
    const auto& banner = parser.get_preamble();

    auto btype = banner.type();
    if (btype == eminem::Field::REAL || btype == eminem::Field::DOUBLE || btype == eminem::Field::INTEGER) {
        Rcpp::NumericVector values(ntotal);

        if (btype == eminem::Field::INTEGER) {
            scan_real([&](eminem::Index r, eminem::Index c, double val) -> void {
                auto& pos = offsets[c];
                row_indices[pos] = r;
                values[pos] = val;
                ++pos;
            });
        } else {
            scan_integer([&](eminem::Index r, eminem::Index c, int val) -> void {
                auto& pos = offsets[c];
                row_indices[pos] = r;
                values[pos] = val;
                ++pos;
            });
        }

        int* iptr = row_indices.begin();
        double* vptr = values.begin();
        const int* ptr = indptr.begin();
        subpar::parallelize_range(threads, NC, [&](int, decltype(NC) start, decltype(NC) length) -> void {
            std::vector<std::pair<int, double> > sortbuffer;
            for (decltype(start) c = start, end = start + length; c < end; ++c) {
                auto pstart = pptr[c], pend = ptr[c + 1];
                if (std::is_sorted(iptr + pstart, iptr + pend)) {
                    continue;
                }
                sortbuffer.clear();
                for (auto p = pstart; p < pend; ++p) {
                    sortbuffer.emplace_back(iptr[p], vptr[p]);
                }
                std::sort(sortbuffer.begin(), sortbuffer.end());
                for (decltype(sortbuffer.size()) i = 0, end = sortbuffer.size(); i < end; ++i) {
                    auto offset = pstart + i;
                    const auto& current = sortbuffer[i];
                    iptr[offset] = current.first;
                    vptr[offset] = current.second;
                }
            }
        });

        return Rcpp::List(
            Rcpp::Named("i") = row_indices, 
            Rcpp::Named("x") = values, 
            Rcpp::Named("p") = indptr
        );

    } else if (btype == eminem::Field::PATTERN) {
        scan_pattern([&](eminem::Index r, eminem::Index c, int val) -> void {
            auto& pos = offsets[c];
            row_indices[pos] = r;
            values[pos] = val;
            ++pos;
        });

        int* iptr = row_indices.begin();
        const int* ptr = indptr.begin();
        subpar::parallelize_range(threads, NC, [&](int, decltype(NC) start, decltype(NC) length) -> void {
            for (decltype(start) c = start, end = start + length; c < end; ++c) {
                auto pstart = pptr[c], pend = ptr[c + 1];
                if (std::is_sorted(iptr + pstart, iptr + pend)) {
                    continue;
                }
                std::sort(iptr + pstart, iptr + pend);
            }
        });

        return Rcpp::List(
            Rcpp::Named("i") = row_indices, 
            Rcpp::Named("p") = indptr
        );

    } else {
        throw std::runtime_error("unknown eminem::Field type");
        return R_NilValue;
    }
}

Rcpp::RObject read_mm_two_pass(const std::string& path, const std::string& class_name, int threads) {
    // First pass, to determine the offsets.
    std::vector<int> nnz_per_col;
    Rcpp::IntegerVector dimensions(2);
    {
        byteme::ParseSomeFileOptions opt;
        opt.threads = threads;
        auto parser = parse_some_file(path.c_str(), opt);
        dimensions[0] = parser.nrow();
        dimensions[1] = parser.ncol();
        nnz_per_col.resize(NC);

        const auto& banner = parser.get_preamble();
        switch (banner.type()) {
            case eminem::Field::REAL:
            case eminem::Field::DOUBLE:
                scan_real([&](eminem::Index, eminem::Index c, double) -> void {
                    ++(offsets[c + 1]);
                });
                break;
            case eminem::Field::INTEGER:
                scan_real([&](eminem::Index, eminem::Index c, int) -> void {
                    ++(offsets[c + 1]);
                });
                break;
            case eminem::Field::PATTERN:
                scan_complex([&](eminem::Index, eminem::Index c) -> void {
                    ++(offsets[c + 1]);
                });
                break;
            default:
                throw std::runtime_error("unknown eminem::Field type");
        }
    }

    // Second pass, to fill the vectors.
    if (class_name == "SVT_SparseMatrix") {
        return Rcpp::List(
            Rcpp::Named("dim") = dimensions,
            Rcpp::Named("contents") = read_mm_two_pass_SVT_SparseMatrix(path, nnz_per_col, threads)
        );
    } else {
        return Rcpp::List(
            Rcpp::Named("dim") = dimensions,
            Rcpp::Named("contents") = read_mm_two_pass_CsparseMatrix(path, nnz_per_col, threads)
        );
    }
}


