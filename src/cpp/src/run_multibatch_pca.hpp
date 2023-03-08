#ifndef SCRANJL_RUN_MULTIBATCH_PCA_HPP
#define SCRANJL_RUN_MULTIBATCH_PCA_HPP

#include "jlcxx/jlcxx.hpp"
#include "ScranMatrix.h"
#include <cstdint>

#ifndef JL_BINDINGS_ONLY
#include "scran/scran.hpp"
#endif

#include "pca_utils.h"

typedef RunPcaResultsBase<scran::MultiBatchPCA::Results> MultiBatchPcaResults;

MultiBatchPcaResults run_multibatch_pca(const ScranMatrix& x, int ndim, bool use_subset, jlcxx::ArrayRef<uint8_t> features, jlcxx::ArrayRef<int32_t> block, bool scale, int nthreads) 
#ifndef JL_BINDINGS_ONLY
{
    scran::MultiBatchPCA pcs;
    pcs
        .set_rank(ndim)
        .set_scale(scale)
        .set_num_threads(nthreads);

    const uint8_t* fptr = (use_subset ? features.data() : NULL);
    auto res = pcs.run(x.ptr.get(), block.data(), fptr);

    return MultiBatchPcaResults(std::move(res));
}
#endif
;

#endif
