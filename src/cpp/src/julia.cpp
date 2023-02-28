#include "jlcxx/jlcxx.hpp"

#include "ScranMatrix.h"
#include "suggest_rna_qc_filters.h"

ScranMatrix initialize_from_memory_int32_int32(jlcxx::ArrayRef<int32_t> x, jlcxx::ArrayRef<int32_t> i, jlcxx::ArrayRef<int32_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix initialize_from_memory_int64_int32(jlcxx::ArrayRef<int64_t> x, jlcxx::ArrayRef<int32_t> i, jlcxx::ArrayRef<int32_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix initialize_from_memory_float32_int32(jlcxx::ArrayRef<float> x, jlcxx::ArrayRef<int32_t> i, jlcxx::ArrayRef<int32_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix initialize_from_memory_float64_int32(jlcxx::ArrayRef<double> x, jlcxx::ArrayRef<int32_t> i, jlcxx::ArrayRef<int32_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix initialize_from_memory_int32_int64(jlcxx::ArrayRef<int32_t> x, jlcxx::ArrayRef<int64_t> i, jlcxx::ArrayRef<int64_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix initialize_from_memory_int64_int64(jlcxx::ArrayRef<int64_t> x, jlcxx::ArrayRef<int64_t> i, jlcxx::ArrayRef<int64_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix initialize_from_memory_float32_int64(jlcxx::ArrayRef<float> x, jlcxx::ArrayRef<int64_t> i, jlcxx::ArrayRef<int64_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix initialize_from_memory_float64_int64(jlcxx::ArrayRef<double> x, jlcxx::ArrayRef<int64_t> i, jlcxx::ArrayRef<int64_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced);

ScranMatrix log_norm_counts(const ScranMatrix& x, 
    bool use_block, 
    jlcxx::ArrayRef<int32_t> block_ids, 
    std::string block_method,
    bool use_size_factors, 
    jlcxx::ArrayRef<double> size_factors, 
    bool center, 
    bool allow_zeros,
    int nthreads);

void per_cell_rna_qc_metrics(const ScranMatrix& x,
    jlcxx::ArrayRef<jl_value_t*> subsets,
    jlcxx::ArrayRef<double> sums,
    jlcxx::ArrayRef<int32_t> detected,
    jlcxx::ArrayRef<jl_value_t*> proportions,
    int nthreads);

SuggestRnaQcFiltersResults suggest_rna_qc_filters(
    jlcxx::ArrayRef<double> sums, 
    jlcxx::ArrayRef<int32_t> detected, 
    jlcxx::ArrayRef<jl_value_t*> proportions, 
    bool use_block, 
    jlcxx::ArrayRef<int32_t> block_ids, 
    jlcxx::ArrayRef<uint8_t> keep, 
    double nmads);

JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
    mod.add_type<ScranMatrix>("ScranMatrix")
       .method("row", &ScranMatrix::row)
       .method("column", &ScranMatrix::column)
       .method("num_rows", &ScranMatrix::nrow)
       .method("num_columns", &ScranMatrix::ncol);

    mod.add_type<SuggestRnaQcFiltersResults>("SuggestRnaQcFiltersResults")
       .method("thresholds_sums", &SuggestRnaQcFiltersResults::thresholds_sums)
       .method("thresholds_detected", &SuggestRnaQcFiltersResults::thresholds_detected)
       .method("thresholds_proportions", &SuggestRnaQcFiltersResults::thresholds_proportions);

    mod.method("initialize_from_memory_int32_int32", &initialize_from_memory_int32_int32);
    mod.method("initialize_from_memory_int64_int32", &initialize_from_memory_int64_int32);
    mod.method("initialize_from_memory_float32_int32", &initialize_from_memory_float32_int32);
    mod.method("initialize_from_memory_float64_int32", &initialize_from_memory_float64_int32);

    mod.method("initialize_from_memory_int32_int64", &initialize_from_memory_int32_int64);
    mod.method("initialize_from_memory_int64_int64", &initialize_from_memory_int64_int64);
    mod.method("initialize_from_memory_float32_int64", &initialize_from_memory_float32_int64);
    mod.method("initialize_from_memory_float64_int64", &initialize_from_memory_float64_int64);

    mod.method("log_norm_counts", &log_norm_counts);

    mod.method("per_cell_rna_qc_metrics", &per_cell_rna_qc_metrics);

    mod.method("suggest_rna_qc_filters", &suggest_rna_qc_filters);
}
