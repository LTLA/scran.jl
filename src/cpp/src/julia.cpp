#include "jlcxx/jlcxx.hpp"

#include "ScranMatrix.h"

ScranMatrix initialize_from_memory_int32_int32(jlcxx::ArrayRef<int32_t> x, jlcxx::ArrayRef<int32_t> i, jlcxx::ArrayRef<int32_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix initialize_from_memory_int64_int32(jlcxx::ArrayRef<int64_t> x, jlcxx::ArrayRef<int32_t> i, jlcxx::ArrayRef<int32_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix initialize_from_memory_float32_int32(jlcxx::ArrayRef<float> x, jlcxx::ArrayRef<int32_t> i, jlcxx::ArrayRef<int32_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix initialize_from_memory_float64_int32(jlcxx::ArrayRef<double> x, jlcxx::ArrayRef<int32_t> i, jlcxx::ArrayRef<int32_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix initialize_from_memory_int32_int64(jlcxx::ArrayRef<int32_t> x, jlcxx::ArrayRef<int64_t> i, jlcxx::ArrayRef<int64_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix initialize_from_memory_int64_int64(jlcxx::ArrayRef<int32_t> x, jlcxx::ArrayRef<int64_t> i, jlcxx::ArrayRef<int64_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix initialize_from_memory_float32_int64(jlcxx::ArrayRef<float> x, jlcxx::ArrayRef<int64_t> i, jlcxx::ArrayRef<int64_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix initialize_from_memory_float64_int64(jlcxx::ArrayRef<double> x, jlcxx::ArrayRef<int64_t> i, jlcxx::ArrayRef<int64_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix log_norm_counts(const ScranMatrix& x, 
    bool use_block, 
    const jlcxx::ArrayRef<int32_t>& block_ids, 
    std::string block_method,
    bool use_size_factors, 
    const jlcxx::ArrayRef<double>& size_factors, 
    bool center, 
    bool allow_zeros,
    int nthreads);

JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
    mod.add_type<ScranMatrix>("ScranMatrix")
       .method("row", &ScranMatrix::row)
       .method("column", &ScranMatrix::column)
       .method("num_rows", &ScranMatrix::nrow)
       .method("num_columns", &ScranMatrix::ncol);

    mod.method("initialize_from_memory_int32_int32", &initialize_from_memory_int32_int32);
    mod.method("initialize_from_memory_int64_int32", &initialize_from_memory_int64_int32);
    mod.method("initialize_from_memory_float32_int32", &initialize_from_memory_float32_int32);
    mod.method("initialize_from_memory_float64_int32", &initialize_from_memory_float64_int32);

    mod.method("initialize_from_memory_int32_int64", &initialize_from_memory_int32_int64);
    mod.method("initialize_from_memory_int64_int64", &initialize_from_memory_int64_int64);
    mod.method("initialize_from_memory_float32_int64", &initialize_from_memory_float32_int64);
    mod.method("initialize_from_memory_float64_int64", &initialize_from_memory_float64_int64);

    mod.method("log_norm_counts", &log_norm_counts);
}
