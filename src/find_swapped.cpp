#include "DropletUtils.h"
#include <cstdint> 

SEXP find_swapped(SEXP _groups, SEXP _reads, SEXP _minfrac) {
    BEGIN_RCPP

    // Checking all of the inputs.
    Rcpp::IntegerVector groups(_groups);
    Rcpp::IntegerVector reads(_reads);
    const int nmolecules=std::accumulate(groups.begin(), groups.end(), 0);
    if (nmolecules!=reads.size()) {
        throw std::runtime_error("length of 'reads' vector should be equal to sum of RLE lengths");
    }

    Rcpp::NumericVector minfrac(_minfrac);
    if (minfrac.size()!=1) { 
        throw std::runtime_error("minimum fraction should be a numeric scalar");
    }
    const double mf=minfrac[0];

    // Setting up the output.
    Rcpp::LogicalVector output(nmolecules, 1);
    auto oIt=output.begin();

    // Iterating across the molecule groups.
    auto rIt=reads.begin();
    for (const auto& g : groups) {
        int topindex=0, topreads=0, totalreads=0;

        for (int i=0; i<g; ++i, ++rIt) {
            if (topreads < *rIt) {
                topreads=*rIt;
                topindex=i;
            }
            totalreads+=*rIt;
        }

        if (double(totalreads)*mf <= double(topreads)) {
            *(oIt+topindex)=0;
            oIt+=g;
        }
    }

    return output;
    END_RCPP;
}

SEXP get_cell_barcodes(SEXP _fname, SEXP _dname, SEXP _barcodelen) {
    BEGIN_RCPP
    Rcpp::StringVector fname(_fname);
    Rcpp::StringVector dname(_dname);
    if (fname.size()!=1) {
        throw std::runtime_error("file name should be a string");
    }
    if (dname.size()!=1) {
        throw std::runtime_error("dataset name should be a string");
    }

    Rcpp::IntegerVector barcodelen(_barcodelen);
    if (barcodelen.size()!=1) {
        throw std::runtime_error("barcode length should be an integer vector");
    }
    const int blen=barcodelen[0];

    // Setting the file input parameters.
    std::string curfile=Rcpp::as<std::string>(fname[0]);
    std::string curdata=Rcpp::as<std::string>(dname[0]);
    H5::H5File h5file(curfile.c_str(), H5F_ACC_RDONLY);
    H5::DataSet h5data = h5file.openDataSet(curdata.c_str());
    
    if (h5data.getTypeClass()!=H5T_INTEGER) {
        throw std::runtime_error("cell barcodes should be encoded as integers");
    }

    // Determining the dimensions.
    H5::DataSpace dataspace = h5data.getSpace();
    if (dataspace.getSimpleExtentNdims()!=1) {
        throw std::runtime_error("cell barcodes should be a one-dimensional array");
    }
    hsize_t dims_out;
    dataspace.getSimpleExtentDims(&dims_out, NULL);

    // Extracting the data.
    H5::DataSpace memspace(1, &dims_out);
    memspace.selectAll();
    dataspace.selectAll();
    std::vector<uint64_t> encoded(dims_out);
    h5data.read(encoded.data(), H5::PredType::NATIVE_INT64, memspace, dataspace);
   
    // Iterating across the output and taking pairs of bits.
    Rcpp::StringVector output(dims_out);
    auto oIt=output.begin();
    std::vector<char> ref(blen+1, '\0');   
    const char* bases="ACGT";

    for (const auto& enc : encoded) {
        for (int pos=0; pos<blen; ++pos) {
            auto leftover=(enc >> 2*pos) & 0x3;
            ref[blen - pos - 1]=bases[leftover];
        }

        (*oIt)=Rcpp::String(ref.data());
        ++oIt;
    }

    return output;    
    END_RCPP
}


