#include "Rcpp.h"

#include "H5Cpp.h"
#include "utils.h"

#include <stdexcept>
#include <cstdint>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

//[[Rcpp::export(rng=false)]]
Rcpp::StringVector get_cell_barcodes(std::string fname, std::string dname, Rcpp::RObject barcodelen) {
    // Setting the file input parameters.
    H5::H5File h5file(fname.c_str(), H5F_ACC_RDONLY);
    H5::DataSet h5data = h5file.openDataSet(dname.c_str());
    
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
    std::vector<std::uint64_t> encoded(dims_out);
    h5data.read(encoded.data(), H5::PredType::NATIVE_UINT64, memspace, dataspace);
   
    // Guessing the barcode length. 
    int blen=0;
    if (barcodelen.isNULL()) {
        if (encoded.size()) { 
            blen=std::ceil(std::log(*std::max_element(encoded.begin(), encoded.end()))/std::log(4));
        }
    } else {
        blen=Rcpp::as<int>(barcodelen);
    }

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
}
