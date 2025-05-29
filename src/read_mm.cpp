#include "eminem/eminem.hpp"
#include "byteme/byteme.hpp"
#include "subpar/subpar.hpp"

#include "Rcpp.h"

#include <vector>
#include <stdexcept>
#include <limits>
#include <type_traits>

template<typename Type_>
void sort_SVT_SparseMatrix_columns(const std::vector<int*>& iptrs, const std::vector<Type_*>& vptrs, const std::vector<int>& num, int threads) {
    auto NC = iptrs.size();
    subpar::parallelize_range(threads, NC, [&](int, decltype(NC) start, decltype(NC) length) -> void {
        std::vector<std::pair<int, Type_> > sortbuffer;
        for (decltype(start) c = start, end = start + length; c < end; ++c) {
            auto iptr = iptrs[c];
            auto n = num[c];
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

Rcpp::RObject read_mm_two_pass_SVT_SparseMatrix(const std::string& path, const std::vector<int>& nnz_per_col, int threads) {
    auto NC = nnz_per_col.size();
    Rcpp::List contents(NC);
    std::vector<int*> iptrs(NC);
    std::vector<int> used(NC);

    eminem::ParseSomeFileOptions opt;
    opt.num_threads = threads;
    auto parser = eminem::parse_some_file(path.c_str(), opt);
    parser.scan_preamble();
    const auto& banner = parser.get_banner();

    if (banner.field == eminem::Field::REAL || banner.field == eminem::Field::DOUBLE) {
        std::vector<double*> vptrs(NC);
        for (decltype(NC) c = 0; c < NC; ++c) {
            Rcpp::NumericVector values(nnz_per_col[c]);
            Rcpp::IntegerVector indices(nnz_per_col[c]);
            iptrs[c] = indices.begin(); // these pointers should still be valid after the std::move as they refer to R-managed allocations.
            vptrs[c] = values.begin();
            contents[c] = Rcpp::List::create(std::move(indices), std::move(values));
        }

        parser.scan_real([&](eminem::Index r, eminem::Index c, double val) -> void {
            auto& pos = used[c - 1];
            iptrs[c - 1][pos] = r - 1;
            vptrs[c - 1][pos] = val;
            ++pos;
        });

        sort_SVT_SparseMatrix_columns(iptrs, vptrs, used, threads);

    } else if (banner.field == eminem::Field::INTEGER) {
        std::vector<int*> vptrs(NC);
        for (decltype(NC) c = 0; c < NC; ++c) {
            Rcpp::IntegerVector values(nnz_per_col[c]);
            Rcpp::IntegerVector indices(nnz_per_col[c]);
            iptrs[c] = indices.begin(); // these pointers should still be valid after the std::move as they refer to R-managed allocations.
            vptrs[c] = values.begin();
            contents[c] = Rcpp::List::create(std::move(indices), std::move(values));
        }

        parser.scan_integer([&](eminem::Index r, eminem::Index c, double val) -> void {
            auto& pos = used[c - 1];
            iptrs[c - 1][pos] = r - 1;
            vptrs[c - 1][pos] = val;
            ++pos;
        });

        sort_SVT_SparseMatrix_columns(iptrs, vptrs, used, threads);

    } else {
        throw std::runtime_error("unsupported eminem::Field type");
    }

    return contents;
}

template<typename Size_>
int safe_add_indptr(int sofar, Size_ val) {
    constexpr auto limiter = std::numeric_limits<int>::max();
    if (static_cast<unsigned>(limiter) < val || static_cast<int>(limiter - val) < sofar) {
        throw std::runtime_error("too many non-zero elements to be stored in a CsparseMatrix");
    }
    return sofar + val;
}

template<typename Output_, typename Number_>
Output_ safe_get_indptr_size(Number_ n) {
    if (static_cast<typename std::make_unsigned<Number_>::type>(n) >= static_cast<typename std::make_unsigned<Output_>::type>(std::numeric_limits<Output_>::max())) {
        throw std::runtime_error("number of columns is too large for allocating indptrs");
    }
    Output_ out = n;
    ++out;
    return out;
}

Rcpp::RObject read_mm_two_pass_CsparseMatrix(const std::string& path, const std::vector<int>& nnz_per_col, int threads) {
    auto NC = nnz_per_col.size();
    std::vector<int> offsets(safe_get_indptr_size<typename std::vector<int>::size_type>(NC));
    for (decltype(NC) c = 0; c < NC; ++c) {
        offsets[c + 1] = safe_add_indptr(offsets[c], nnz_per_col[c]);
    }

    Rcpp::IntegerVector indptr(offsets.begin(), offsets.end());
    auto ntotal = indptr[NC];
    Rcpp::IntegerVector row_indices(ntotal);

    eminem::ParseSomeFileOptions opt;
    opt.num_threads = threads;
    auto parser = eminem::parse_some_file(path.c_str(), opt);
    parser.scan_preamble();
    const auto& banner = parser.get_banner();

    if (banner.field == eminem::Field::REAL || banner.field == eminem::Field::DOUBLE || banner.field == eminem::Field::INTEGER) {
        Rcpp::NumericVector values(ntotal);

        if (banner.field == eminem::Field::INTEGER) {
            parser.scan_real([&](eminem::Index r, eminem::Index c, double val) -> void {
                auto& pos = offsets[c];
                row_indices[pos] = r;
                values[pos] = val;
                ++pos;
            });
        } else {
            parser.scan_integer([&](eminem::Index r, eminem::Index c, int val) -> void {
                auto& pos = offsets[c];
                row_indices[pos] = r;
                values[pos] = val;
                ++pos;
            });
        }

        int* iptr = row_indices.begin();
        double* vptr = values.begin();
        const int* pptr = indptr.begin();
        subpar::parallelize_range(threads, NC, [&](int, decltype(NC) start, decltype(NC) length) -> void {
            std::vector<std::pair<int, double> > sortbuffer;
            for (decltype(start) c = start, end = start + length; c < end; ++c) {
                auto pstart = pptr[c], pend = pptr[c + 1];
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

        return Rcpp::List::create(
            Rcpp::Named("i") = row_indices, 
            Rcpp::Named("x") = values, 
            Rcpp::Named("p") = indptr
        );

    } else {
        throw std::runtime_error("unsupported eminem::Field type");
        return R_NilValue;
    }
}

template<typename Size_>
int safe_cast_dim(Size_ val) {
    constexpr auto limiter = std::numeric_limits<int>::max();
    if (static_cast<unsigned>(limiter) < val) {
        throw std::runtime_error("dimension extent is too large to be stored as an integer");
    }
    return val;
}

Rcpp::RObject read_mm_two_pass(const std::string& path, const std::string& class_name, int threads) {
    // First pass, to determine the size of each column for preallocation.
    std::vector<int> nnz_per_col;
    Rcpp::IntegerVector dimensions(2);
    eminem::ParseSomeFileOptions opt;
    opt.num_threads = threads;
    auto parser = eminem::parse_some_file(path.c_str(), opt);
    parser.scan_preamble();

    dimensions[0] = safe_cast_dim(parser.get_nrows());
    auto NC = safe_cast_dim(parser.get_ncols());
    dimensions[1] = NC;

    nnz_per_col.resize(NC);
    const auto& banner = parser.get_banner();
    switch (banner.field) {
        case eminem::Field::REAL: case eminem::Field::DOUBLE:
            // Don't bother checking for overflow, as we already did that for the dimension
            // extents and eminem will automatically check that indices lie within range.
            // Note that indices are 1-based. 
            parser.scan_real([&](eminem::Index, eminem::Index c, double) -> void {
                ++(nnz_per_col[c - 1]);
            });
            break;
        case eminem::Field::INTEGER:
            parser.scan_real([&](eminem::Index, eminem::Index c, int) -> void {
                ++(nnz_per_col[c - 1]);
            });
            break;
        default:
            throw std::runtime_error("unsupported eminem::Field type");
    }

    // Second pass, to fill the vectors.
    if (class_name == "SVT_SparseMatrix") {
        return Rcpp::List::create(
            Rcpp::Named("dim") = dimensions,
            Rcpp::Named("contents") = read_mm_two_pass_SVT_SparseMatrix(path, nnz_per_col, threads)
        );
    } else {
        return Rcpp::List::create(
            Rcpp::Named("dim") = dimensions,
            Rcpp::Named("contents") = read_mm_two_pass_CsparseMatrix(path, nnz_per_col, threads)
        );
    }
}

template<typename Rclass_, typename Type_>
Rcpp::RObject format_one_pass_output(std::vector<std::pair<std::vector<int>, std::vector<Type_> > >& contents, const std::string& class_name, int threads) {
    auto NC = contents.size();
    subpar::parallelize_range(threads, NC, [&](int, decltype(NC) start, decltype(NC) length) -> void {
        std::vector<std::pair<int, Type_> > sortbuffer;
        for (decltype(start) c = start, end = start + length; c < end; ++c) {
            auto& idxs = contents[c].first;
            if (std::is_sorted(idxs.begin(), idxs.end())) {
                continue;
            }
            auto& vals = contents[c].second;
            auto n = idxs.size();
            sortbuffer.clear();
            for (decltype(n) i = 0; i < n; ++i) {
                sortbuffer.emplace_back(idxs[i], vals[i]);
            }
            std::sort(sortbuffer.begin(), sortbuffer.end());
            for (decltype(n) i = 0; i < n; ++i) {
                const auto& current = sortbuffer[i];
                idxs[i] = current.first;
                vals[i] = current.second;
            }
        }
    });

    if (class_name == "SVT_SparseMatrix") {
        Rcpp::List output(NC);
        for (decltype(NC) c = 0; c < NC; ++c) {
            const auto& pair = contents[c];
            output[c] = Rcpp::List::create(
                Rcpp::IntegerVector(pair.first.begin(), pair.first.end()),
                Rclass_(pair.second.begin(), pair.second.end())
            );
        }
        return output;

    } else {
        Rcpp::IntegerVector indptr(safe_get_indptr_size<R_xlen_t>(NC));
        for (decltype(NC) c = 0; c < NC; ++c) {
            indptr[c + 1] = safe_add_indptr(indptr[c], contents[c].first.size());
        }

        auto total_nnz = indptr[NC];
        Rcpp::IntegerVector indices(total_nnz);
        Rcpp::NumericVector values(total_nnz); // it's going to be a dgCMatrix anyway, so we might as well save it as a numeric vector.
        decltype(total_nnz) sofar = 0; 
        for (decltype(NC) c = 0; c < NC; ++c) {
            const auto& pair = contents[c];
            std::copy(pair.first.begin(), pair.first.end(), indices.begin() + sofar);
            std::copy(pair.second.begin(), pair.second.end(), values.begin() + sofar);
            sofar += pair.first.size(); 
        }

        return Rcpp::List::create(
            Rcpp::Named("i") = indices,
            Rcpp::Named("x") = values,
            Rcpp::Named("p") = indptr
        );
    }
}

Rcpp::RObject read_mm_one_pass(const std::string& path, const std::string& class_name, int threads) {
    eminem::ParseSomeFileOptions opt;
    opt.num_threads = threads;
    auto parser = eminem::parse_some_file(path.c_str(), opt);
    parser.scan_preamble();

    Rcpp::IntegerVector dimensions(2);
    dimensions[0] = safe_cast_dim(parser.get_nrows());
    auto NC = safe_cast_dim(parser.get_ncols());
    dimensions[1] = NC;

    const auto& banner = parser.get_banner();

    if (banner.field == eminem::Field::REAL || banner.field == eminem::Field::DOUBLE) {
        std::vector<std::pair<std::vector<int>, std::vector<double> > > contents(NC);
        parser.scan_real([&](eminem::Index r, eminem::Index c, double val) -> void {
            contents[c - 1].first.push_back(r - 1);
            contents[c - 1].second.push_back(val);
        });
        return Rcpp::List::create(
            Rcpp::Named("dim") = dimensions,
            Rcpp::Named("contents") = format_one_pass_output<Rcpp::NumericVector>(contents, class_name, threads)
        );

    } else if (banner.field == eminem::Field::INTEGER) {
        std::vector<std::pair<std::vector<int>, std::vector<int> > > contents(NC);
        parser.scan_real([&](eminem::Index r, eminem::Index c, int val) -> void {
            contents[c - 1].first.push_back(r - 1);
            contents[c - 1].second.push_back(val);
        });
        return Rcpp::List::create(
            Rcpp::Named("dim") = dimensions,
            Rcpp::Named("contents") = format_one_pass_output<Rcpp::IntegerVector>(contents, class_name, threads)
        );

    } else {
        throw std::runtime_error("unsupported eminem::Field type");
        return R_NilValue;
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::RObject read_mm(const std::string& path, bool two_pass, const std::string& class_name, int threads) {
    if (two_pass) {
        return read_mm_two_pass(path, class_name, threads);
    } else {
        return read_mm_one_pass(path, class_name, threads);
    }
}
