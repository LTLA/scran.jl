#include "jlcxx/jlcxx.hpp"

#include "ScranMatrix.h"
#include "scran/scran.hpp"

#include <vector>
#include <cstdint>
#include <limits>
#include <type_traits>

ScranMatrix log_norm_counts(const ScranMatrix& x, 
    bool use_block, 
    jlcxx::ArrayRef<int32_t> block_ids, 
    std::string block_method,
    bool use_size_factors, 
    jlcxx::ArrayRef<double> size_factors, 
    bool center, 
    bool allow_zeros,
    int nthreads)
{
    scran::LogNormCounts logger;
    logger
        .set_center(center)
        .set_handle_zeros(allow_zeros)
        .set_num_threads(nthreads);

    if (block_method == "lowest") {
        logger.set_block_mode(scran::CenterSizeFactors::LOWEST);
    } else if (block_method == "perblock") {
        logger.set_block_mode(scran::CenterSizeFactors::PER_BLOCK);
    } else {
        throw std::runtime_error("unknown choice '" + block_method + "' for the blocking method");
    }

    const int32_t* bptr = (use_block ? block_ids.data() : NULL);
    if (use_size_factors) {
        return ScranMatrix(logger.run_blocked(x.ptr, std::vector<double>(size_factors.begin(), size_factors.end()), bptr));
    } else {
        return ScranMatrix(logger.run_blocked(x.ptr, bptr));
    }
}
