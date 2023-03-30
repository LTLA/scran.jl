#ifndef SCRANJL_GROUPED_SIZE_FACTORS_HPP
#define SCRANJL_GROUPED_SIZE_FACTORS_HPP

#include "jlcxx/jlcxx.hpp"
#include "ScranMatrix.h"
#include <cstdint>

#ifndef JL_BINDINGS_ONLY
#include "scran/scran.hpp"
#include "scran/normalization/GroupedSizeFactors.hpp"
#endif

void grouped_size_factors(
    const ScranMatrix& x, 
    jlcxx::ArrayRef<int32_t> cluster_ids, 
    jlcxx::ArrayRef<double> output, 
    bool center, 
    double prior_count, 
    int reference, 
    int nthreads) 
#ifndef JL_BINDINGS_ONLY
{
    scran::GroupedSizeFactors runner;
    runner.set_center(center).set_prior_count(prior_count).set_num_threads(nthreads);

    const int32_t* cptr = cluster_ids.data();
    double* optr = output.data();
    if (reference >= 0) {
        runner.run(x.ptr.get(), cptr, reference, optr);
    } else {
        runner.run(x.ptr.get(), cptr, optr);
    }

    return;
}
#endif
;

#endif
