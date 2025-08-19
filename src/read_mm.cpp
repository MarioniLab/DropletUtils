#include "eminem/eminem.hpp"
#include "byteme/byteme.hpp"
#include "subpar/subpar.hpp"
#include "sanisizer/sanisizer.hpp"

#include "Rcpp.h"

#include <vector>
#include <stdexcept>
#include <limits>
#include <type_traits>
#include <algorithm>

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
    auto contents = sanisizer::create<Rcpp::List>(NC);
    auto iptrs = sanisizer::create<std::vector<int*> >(NC);
    auto used = sanisizer::create<std::vector<int> >(NC);

    auto parser = eminem::parse_some_file<int>(path.c_str(), [&]{
        eminem::ParseSomeFileOptions opt;
        opt.num_threads = threads;
        return opt;
    }());
    parser.scan_preamble();
    const auto& banner = parser.get_banner();
    std::string out_type;

    if (banner.field == eminem::Field::REAL || banner.field == eminem::Field::DOUBLE) {
        auto vptrs = sanisizer::create<std::vector<double*> >(NC);
        for (decltype(NC) c = 0; c < NC; ++c) {
            auto values = sanisizer::create<Rcpp::NumericVector>(nnz_per_col[c]);
            auto indices = sanisizer::create<Rcpp::IntegerVector>(nnz_per_col[c]);
            iptrs[c] = indices.begin(); // these pointers should still be valid after the std::move as they refer to R-managed allocations.
            vptrs[c] = values.begin();
            contents[c] = Rcpp::List::create(std::move(values), std::move(indices));
        }

        parser.scan_real([&](int r, int c, double val) -> void {
            auto& pos = used[c - 1];
            iptrs[c - 1][pos] = r - 1;
            vptrs[c - 1][pos] = val;
            ++pos;
        });

        sort_SVT_SparseMatrix_columns(iptrs, vptrs, used, threads);
        out_type = "double";

    } else if (banner.field == eminem::Field::INTEGER) {
        auto vptrs = sanisizer::create<std::vector<int*> >(NC);
        for (decltype(NC) c = 0; c < NC; ++c) {
            auto values = sanisizer::create<Rcpp::IntegerVector>(nnz_per_col[c]);
            auto indices = sanisizer::create<Rcpp::IntegerVector>(nnz_per_col[c]);
            iptrs[c] = indices.begin(); // these pointers should still be valid after the std::move as they refer to R-managed allocations.
            vptrs[c] = values.begin();
            contents[c] = Rcpp::List::create(std::move(values), std::move(indices));
        }

        parser.scan_integer([&](int r, int c, double val) -> void {
            auto& pos = used[c - 1];
            iptrs[c - 1][pos] = r - 1;
            vptrs[c - 1][pos] = val;
            ++pos;
        });

        sort_SVT_SparseMatrix_columns(iptrs, vptrs, used, threads);
        out_type = "integer";

    } else {
        throw std::runtime_error("unsupported eminem::Field type");
    }

    return Rcpp::List::create(
        Rcpp::Named("list") = contents,
        Rcpp::Named("type") = out_type 
    );
}

Rcpp::RObject read_mm_two_pass_CsparseMatrix(const std::string& path, const std::vector<int>& nnz_per_col, int threads) {
    auto NC = nnz_per_col.size();
    std::vector<int> offsets(sanisizer::sum<typename std::vector<int>::size_type>(NC, 1));
    for (decltype(NC) c = 0; c < NC; ++c) {
        // Increment is safe from overflow as 'c + 1 <= NC'.
        offsets[c + 1] = sanisizer::sum<int>(offsets[c], nnz_per_col[c]);
    }

    Rcpp::IntegerVector indptr(offsets.begin(), offsets.end());
    auto ntotal = indptr[NC];
    auto row_indices = sanisizer::create<Rcpp::IntegerVector>(ntotal);

    auto parser = eminem::parse_some_file<int>(path.c_str(), [&]{
        eminem::ParseSomeFileOptions opt;
        opt.num_threads = threads;
        return opt;
    }());
    parser.scan_preamble();
    const auto& banner = parser.get_banner();

    if (banner.field == eminem::Field::REAL || banner.field == eminem::Field::DOUBLE || banner.field == eminem::Field::INTEGER) {
        auto values = sanisizer::create<Rcpp::NumericVector>(ntotal);

        if (banner.field == eminem::Field::INTEGER) {
            parser.scan_integer([&](int r, int c, int val) -> void {
                auto& pos = offsets[c - 1];
                row_indices[pos] = r - 1;
                values[pos] = val;
                ++pos;
            });
        } else {
            parser.scan_real([&](int r, int c, double val) -> void {
                auto& pos = offsets[c - 1];
                row_indices[pos] = r - 1;
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
                auto pstart = pptr[c], pend = pptr[c + 1]; // increment won't overflow as 'c + 1 <= end'.
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

Rcpp::RObject read_mm_two_pass(const std::string& path, const std::string& class_name, int threads) {
    // First pass, to determine the size of each column for preallocation.
    auto parser = eminem::parse_some_file<int>(path.c_str(), [&]{
        eminem::ParseSomeFileOptions opt;
        opt.num_threads = threads;
        return opt;
    }());
    parser.scan_preamble();

    Rcpp::IntegerVector dimensions(2);
    dimensions[0] = parser.get_nrows();
    auto NC = parser.get_ncols();
    dimensions[1] = NC;

    auto nnz_per_col = sanisizer::create<std::vector<int> >(NC);
    const auto& banner = parser.get_banner();
    switch (banner.field) {
        case eminem::Field::REAL: case eminem::Field::DOUBLE:
            parser.scan_real([&](int, int c, double) -> void {
                sanisizer::sum<int>(nnz_per_col[c - 1], 1);
            });
            break;
        case eminem::Field::INTEGER:
            parser.scan_real([&](int, int c, int) -> void {
                sanisizer::sum<int>(nnz_per_col[c - 1], 1);
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
        auto output = sanisizer::create<Rcpp::List>(NC);
        for (decltype(NC) c = 0; c < NC; ++c) {
            const auto& pair = contents[c];
            output[c] = Rcpp::List::create(
                Rclass_(pair.second.begin(), pair.second.end()),
                Rcpp::IntegerVector(pair.first.begin(), pair.first.end())
            );
        }
        return Rcpp::List::create(
            Rcpp::Named("list") = output,
            Rcpp::Named("type") = []{
                if constexpr(std::is_same<Type_, int>::value) {
                    return std::string("integer");
                } else {
                    return std::string("double");
                }
            }()
        );

    } else {
        Rcpp::IntegerVector indptr(sanisizer::sum<decltype(std::declval<Rcpp::IntegerVector>().size())>(NC, 1));
        for (decltype(NC) c = 0; c < NC; ++c) {
            // Increment is safe from overflow as 'c + 1 <= NC'.
            indptr[c + 1] = sanisizer::sum<int>(indptr[c], contents[c].first.size());
        }

        auto total_nnz = indptr[NC];
        auto indices = sanisizer::create<Rcpp::IntegerVector>(total_nnz);
        auto values = sanisizer::create<Rcpp::NumericVector>(total_nnz); // it's going to be a dgCMatrix anyway, so we might as well save it as a numeric vector.
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
    auto parser = eminem::parse_some_file<int>(path.c_str(), [&]{
        eminem::ParseSomeFileOptions opt;
        opt.num_threads = threads;
        return opt;
    }());
    parser.scan_preamble();

    Rcpp::IntegerVector dimensions(2);
    dimensions[0] = parser.get_nrows();
    auto NC = parser.get_ncols();
    dimensions[1] = NC;

    const auto& banner = parser.get_banner();

    if (banner.field == eminem::Field::REAL || banner.field == eminem::Field::DOUBLE) {
        auto contents = sanisizer::create<std::vector<std::pair<std::vector<int>, std::vector<double> > > >(NC);
        parser.scan_real([&](int r, int c, double val) -> void {
            contents[c - 1].first.push_back(r - 1);
            contents[c - 1].second.push_back(val);
        });
        return Rcpp::List::create(
            Rcpp::Named("dim") = dimensions,
            Rcpp::Named("contents") = format_one_pass_output<Rcpp::NumericVector>(contents, class_name, threads)
        );

    } else if (banner.field == eminem::Field::INTEGER) {
        auto contents = sanisizer::create<std::vector<std::pair<std::vector<int>, std::vector<int> > > >(NC);
        parser.scan_real([&](int r, int c, int val) -> void {
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
