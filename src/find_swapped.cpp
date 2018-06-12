#include "DropletUtils.h"

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

SEXP find_swapped_ultra(SEXP cells, SEXP genes, SEXP umis, SEXP reads, SEXP minfrac) {
    BEGIN_RCPP

    auto Cells=process_list<Rcpp::StringVector>(cells);
    auto Genes=process_list<Rcpp::IntegerVector>(genes);
    auto Umis=process_list<Rcpp::StringVector>(umis);
    auto Reads=process_list<Rcpp::IntegerVector>(reads);

    const double mf=check_numeric_scalar(minfrac, "minimum fraction");
    compare_lists(Cells, Genes);
    compare_lists(Cells, Umis);
    compare_lists(Cells, Reads);

    const size_t nsamples=Cells.size();
    size_t nmolecules=0;
    for (size_t i=0; i<nsamples; ++i) {
        nmolecules+=Cells[i].size();
    } 

    typedef std::pair<size_t, size_t> molecule;
    std::vector<molecule> ordering(nmolecules);
    auto ordIt=ordering.begin();
    for (size_t i=0; i<nsamples; ++i) {
        const size_t cur_nmol=Cells[i].size();
        for (size_t j=0; j<cur_nmol; ++j) {
            ordIt->first = i;
            ordIt->second = j;
            ++ordIt;
        }
    }
    
    // Sorting the indices based on the list values.
    std::sort(ordering.begin(), ordering.end(), [&](const molecule& left, const molecule& right) {
        const auto& left_gene=Genes[left.first][left.second];
        const auto& right_gene=Genes[right.first][right.second];
        if (left_gene < right_gene) {
            return true;
        } else if (left_gene > right_gene) {
            return false;
        }

        const auto& left_umi=Umis[left.first][left.second];
        const auto& right_umi=Umis[right.first][right.second];
        if (left_umi < right_umi) {
            return true;
        } else if (left_umi > right_umi) {
            return false;
        }

        const auto& left_cell=Cells[left.first][left.second];
        const auto& right_cell=Cells[right.first][right.second];
        return left_cell < right_cell;
    });
    
    // Setting up the output, indicating which values to keep from each sample.
    std::vector<Rcpp::LogicalVector> output(nsamples);
    for (size_t i=0; i<nsamples; ++i) {
        output[i]=Rcpp::LogicalVector(Cells[i].size());
    }

    // Iterating across runs of the same UMI/gene/cell combination.
    auto ostart=ordering.begin(), oend=ordering.begin();
    while (ostart!=ordering.end()) {
        int max_nread=Reads[ostart->first][ostart->second];
        int total_nreads=max_nread;
        auto best_mol=ostart;

        while (oend!=ordering.end()
                && Genes[ostart->first][ostart->second]==Genes[oend->first][oend->second] 
                && Umis[ostart->first][ostart->second]==Umis[oend->first][oend->second]
                && Cells[ostart->first][ostart->second]==Cells[oend->first][oend->second]) {
        
            const int current_nread=Reads[oend->first][oend->second];
            if (current_nread > max_nread) {
                max_nread=current_nread;
                best_mol=oend;
            }

            total_nreads += current_nread;
            ++oend;
        }

        if (double(max_nread)/total_nreads >= mf) {
            output[best_mol->first][best_mol->second]=1;
        }
        ostart=oend;
    }

    // Creating the output object.
    Rcpp::List outlist(nsamples);
    for (size_t i=0; i<nsamples; ++i) {
        outlist[i]=output[i];
    }
    return outlist;
    END_RCPP
}


SEXP find_swapped(SEXP _groups, SEXP _reads, SEXP _minfrac) {
    BEGIN_RCPP

    // Checking all of the inputs.
    Rcpp::IntegerVector groups(_groups);
    Rcpp::IntegerVector reads(_reads);
    const int nmolecules=std::accumulate(groups.begin(), groups.end(), 0);
    if (nmolecules!=reads.size()) {
        throw std::runtime_error("length of 'reads' vector should be equal to sum of RLE lengths");
    }
    const double mf=check_numeric_scalar(_minfrac, "minimum fraction");

    // Setting up the output.
    Rcpp::LogicalVector output(nmolecules, 1);
    auto oIt=output.begin();

    // Iterating across the molecule groups.
    auto rIt=reads.begin();
    for (const auto& g : groups) {
        int topindex=0, topreads=*rIt, totalreads=topreads;
        ++rIt;

        for (int i=1; i<g; ++i, ++rIt) {
            if (topreads < *rIt) {
                topreads=*rIt;
                topindex=i;
            }
            totalreads+=*rIt;
        }
        
        if (double(totalreads)*mf <= double(topreads)) {
            *(oIt+topindex)=0;
        }
        oIt+=g;
    }

    return output;
    END_RCPP;
}
