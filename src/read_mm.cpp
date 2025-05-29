#include "eminem/eminem.hpp"
#include "byteme/byteme.hpp"
#include "subpar/subpar.hpp"

#include "Rcpp.h"

#include <vector>
#include <stdexcept>
#include <limits>

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
            auto& pos = used[c];
            iptrs[c][pos] = r;
            vptrs[c][pos] = val;
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

Rcpp::RObject read_mm_two_pass_CsparseMatrix(const std::string& path, const std::vector<int>& nnz_per_col, int threads) {
    auto NC = nnz_per_col.size();
    std::vector<int> offsets(NC + 1);
    constexpr auto limiter = std::numeric_limits<int>::max();
    for (decltype(NC) c = 0; c < NC; ++c) {
        auto curnnz = nnz_per_col[c];
        auto& lastoff = offsets[c];
        if (limiter - curnnz < lastoff) {
            throw std::runtime_error("too many non-zero elements to be stored in a CsparseMatrix");
        }
        offsets[c + 1] = lastoff + curnnz;
    }

    Rcpp::IntegerVector indptr(offsets.begin(), offsets.end());
    auto ntotal = indptr[NC + 1];
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

    } else if (banner.field == eminem::Field::PATTERN) {
        parser.scan_pattern([&](eminem::Index r, eminem::Index c, bool) -> void {
            auto& pos = offsets[c];
            row_indices[pos] = r;
            ++pos;
        });

        int* iptr = row_indices.begin();
        const int* pptr = indptr.begin();
        subpar::parallelize_range(threads, NC, [&](int, decltype(NC) start, decltype(NC) length) -> void {
            for (decltype(start) c = start, end = start + length; c < end; ++c) {
                auto pstart = pptr[c], pend = pptr[c + 1];
                if (std::is_sorted(iptr + pstart, iptr + pend)) {
                    continue;
                }
                std::sort(iptr + pstart, iptr + pend);
            }
        });

        return Rcpp::List::create(
            Rcpp::Named("i") = row_indices, 
            Rcpp::Named("p") = indptr
        );

    } else {
        throw std::runtime_error("unknown eminem::Field type");
        return R_NilValue;
    }
}

Rcpp::RObject read_mm_two_pass(const std::string& path, const std::string& class_name, int threads) {
    // First pass, to determine the size of each column for preallocation.
    std::vector<int> nnz_per_col;
    Rcpp::IntegerVector dimensions(2);
    {
        eminem::ParseSomeFileOptions opt;
        opt.num_threads = threads;
        auto parser = eminem::parse_some_file(path.c_str(), opt);
        parser.scan_preamble();

        dimensions[0] = parser.get_nrows();
        auto NC = parser.get_ncols();
        dimensions[1] = NC;

        nnz_per_col.resize(NC);
        constexpr auto limiter = std::numeric_limits<int>::max();
        auto check_overflow = [&](int count) {
            if (count == limiter) {
                throw std::runtime_error("integer overflow on the number of elements in a column");
            }
        };

        const auto& banner = parser.get_banner();
        switch (banner.field) {
            case eminem::Field::REAL:
            case eminem::Field::DOUBLE:
                parser.scan_real([&](eminem::Index, eminem::Index c, double) -> void {
                    auto& off = nnz_per_col[c + 1];
                    check_overflow(off);
                    ++off;
                });
                break;
            case eminem::Field::INTEGER:
                parser.scan_real([&](eminem::Index, eminem::Index c, int) -> void {
                    auto& off = nnz_per_col[c + 1];
                    check_overflow(off);
                    ++off;
                });
                break;
            case eminem::Field::PATTERN:
                parser.scan_pattern([&](eminem::Index, eminem::Index c, bool) -> void {
                    auto& off = nnz_per_col[c + 1];
                    check_overflow(off);
                    ++off;
                });
                break;
            default:
                throw std::runtime_error("unknown eminem::Field type");
        }
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
        Rcpp::IntegerVector indptr(NC + 1);
        constexpr auto limiter = std::numeric_limits<int>::max();
        for (decltype(NC) c = 0; c < NC; ++c) {
            auto curnnz = contents[c].first.size();
            auto& lastoff = indptr[c];
            if (limiter - curnnz < lastoff) {
                throw std::runtime_error("too many non-zero elements to be stored in a CsparseMatrix");
            }
            indptr[c + 1] = lastoff + curnnz;
        }

        auto total_nnz = indptr[NC + 1];
        Rcpp::IntegerVector indices(total_nnz);
        Rclass_ values(total_nnz);
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
    dimensions[0] = parser.get_nrows();
    auto NC = parser.get_ncols();
    dimensions[1] = NC;

    const auto& banner = parser.get_banner();

    if (banner.field == eminem::Field::REAL || banner.field == eminem::Field::DOUBLE) {
        std::vector<std::pair<std::vector<int>, std::vector<double> > > contents(NC);
        parser.scan_real([&](eminem::Index r, eminem::Index c, double val) -> void {
            contents[c].first.push_back(r);
            contents[c].second.push_back(val);
        });
        return Rcpp::List::create(
            Rcpp::Named("dim") = dimensions,
            Rcpp::Named("contents") = format_one_pass_output<Rcpp::NumericVector>(contents, class_name, threads)
        );

    } else if (banner.field == eminem::Field::INTEGER) {
        std::vector<std::pair<std::vector<int>, std::vector<int> > > contents(NC);
        parser.scan_real([&](eminem::Index r, eminem::Index c, double val) -> void {
            contents[c].first.push_back(r);
            contents[c].second.push_back(val);
        });
        return Rcpp::List::create(
            Rcpp::Named("dim") = dimensions,
            Rcpp::Named("contents") = format_one_pass_output<Rcpp::IntegerVector>(contents, class_name, threads)
        );

    } else if (banner.field == eminem::Field::PATTERN) {
        std::vector<std::vector<int> > contents(NC);
        parser.scan_real([&](eminem::Index r, eminem::Index c, bool) -> void {
            contents[c].push_back(r);
        });

        subpar::parallelize_range(threads, NC, [&](int, decltype(NC) start, decltype(NC) length) -> void {
            for (decltype(start) c = start, end = start + length; c < end; ++c) {
                auto& current = contents[c];
                if (std::is_sorted(current.begin(), current.end())) {
                    continue;
                }
                std::sort(current.begin(), current.end());
            }
        });

        if (class_name == "SVT_SparseMatrix") {
            Rcpp::List output(NC);
            for (decltype(NC) c = 0; c < NC; ++c) {
                const auto& pair = contents[c];
                output[c] = Rcpp::List::create(
                    Rcpp::IntegerVector(pair.begin(), pair.end()),
                    R_NilValue
                );
            }
            return Rcpp::List::create(
                Rcpp::Named("dim") = dimensions,
                Rcpp::Named("contents") = output
            );

        } else {
            Rcpp::IntegerVector indptr(NC + 1);
            constexpr auto limiter = std::numeric_limits<int>::max();
            for (decltype(NC) c = 0; c < NC; ++c) {
                auto curnnz = contents[c].size();
                auto& lastoff = indptr[c];
                if (limiter - curnnz < lastoff) {
                    throw std::runtime_error("more non-zero elements than are supported in a CsparseMatrix");
                }
                indptr[c + 1] = lastoff + curnnz;
            }

            auto total_nnz = indptr[NC + 1];
            Rcpp::IntegerVector indices(total_nnz);
            decltype(total_nnz) sofar = 0; 
            for (decltype(NC) c = 0; c < NC; ++c) {
                const auto& idxs = contents[c];
                std::copy(idxs.begin(), idxs.end(), indices.begin() + sofar);
                sofar += idxs.size(); 
            }

            return Rcpp::List::create(
                Rcpp::Named("dim") = dimensions,
                Rcpp::Named("contents") = Rcpp::List::create(
                    Rcpp::Named("i") = indices,
                    Rcpp::Named("p") = indptr
                )
            );
        }

    } else {
        throw std::runtime_error("unknown eminem::Field type");
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
