#ifndef SCRANJL_RUN_PCA_HPP
#define SCRANJL_RUN_PCA_HPP

#include "jlcxx/jlcxx.hpp"
#include "ScranMatrix.h"
#include <cstdint>

#ifndef JL_BINDINGS_ONLY
#include "scran/scran.hpp"
#endif

struct RunPcaResults {
    RunPcaResults(scran::RunPCA::Results res) : results(std::move(res)) {}

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

    scran::RunPCA::Results results;
};

RunPcaResults run_pca(const ScranMatrix& x, int ndim, bool use_subset, jlcxx::ArrayRef<uint8_t> features, bool scale, int nthreads) 
#ifndef JL_BINDINGS_ONLY
{
    scran::RunPCA pcs;
    pcs
        .set_rank(ndim)
        .set_scale(scale)
        .set_num_threads(nthreads);

    const uint8_t* fptr = (use_subset ? features.data() : NULL);
    auto res = pcs.run(x.ptr.get(), fptr);

    return RunPcaResults(std::move(res));
}
#endif
;

#endif

