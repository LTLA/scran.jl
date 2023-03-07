#ifndef SCRANJL_PCA_UTILS_H
#define SCRANJL_PCA_UTILS_H

#include "jlcxx/jlcxx.hpp"

template<class Results>
struct RunPcaResultsBase {
    RunPcaResultsBase(Results res) : results(std::move(res)) {}

    int num_dims() {
        return results.pcs.rows();
    }

    int num_obs() {
        return results.pcs.cols();
    }

    jlcxx::ArrayRef<double> principal_components() {
        return jlcxx::ArrayRef<double>(results.pcs.data(), results.pcs.rows() * results.pcs.cols());
    }

    int num_genes() {
        return results.rotation.rows();
    }

    jlcxx::ArrayRef<double> rotation() {
        return jlcxx::ArrayRef<double>(results.rotation.data(), results.rotation.rows() * results.rotation.cols());
    }

    jlcxx::ArrayRef<double> variance_explained() {
        return jlcxx::ArrayRef<double>(results.variance_explained.data(), results.variance_explained.size());
    }

    double total_variance() {
        return results.total_variance;
    }

    Results results;
};

#endif
