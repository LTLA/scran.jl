#include "jlcxx/jlcxx.hpp"

#define JL_BINDINGS_ONLY
#include "initialize_from_memory.hpp"
#include "log_norm_counts.hpp"
#include "per_cell_rna_qc_metrics.hpp"
#include "suggest_rna_qc_filters.hpp"
#include "per_cell_adt_qc_metrics.hpp"
#include "suggest_adt_qc_filters.hpp"
#include "per_cell_crispr_qc_metrics.hpp"
#include "suggest_crispr_qc_filters.hpp"
#include "filter_cells.hpp"

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

    mod.add_type<SuggestAdtQcFiltersResults>("SuggestAdtQcFiltersResults")
       .method("thresholds_detected", &SuggestAdtQcFiltersResults::thresholds_detected)
       .method("thresholds_totals", &SuggestAdtQcFiltersResults::thresholds_totals);

    mod.add_type<SuggestCrisprQcFiltersResults>("SuggestCrisprQcFiltersResults")
       .method("thresholds_max_count", &SuggestCrisprQcFiltersResults::thresholds_max_count);

    // initialize_from_memory.hpp
    mod.method("initialize_from_memory_int32_int32", &initialize_from_memory_int32_int32);
    mod.method("initialize_from_memory_int64_int32", &initialize_from_memory_int64_int32);
    mod.method("initialize_from_memory_float32_int32", &initialize_from_memory_float32_int32);
    mod.method("initialize_from_memory_float64_int32", &initialize_from_memory_float64_int32);

    mod.method("initialize_from_memory_int32_int64", &initialize_from_memory_int32_int64);
    mod.method("initialize_from_memory_int64_int64", &initialize_from_memory_int64_int64);
    mod.method("initialize_from_memory_float32_int64", &initialize_from_memory_float32_int64);
    mod.method("initialize_from_memory_float64_int64", &initialize_from_memory_float64_int64);

    // log_norm_counts.hpp
    mod.method("log_norm_counts", &log_norm_counts);

    // per_cell_rna_qc_metrics.hpp
    mod.method("per_cell_rna_qc_metrics", &per_cell_rna_qc_metrics);

    // suggest_rna_qc_filters.hpp
    mod.method("suggest_rna_qc_filters", &suggest_rna_qc_filters);

    // per_cell_adt_qc_metrics.hpp
    mod.method("per_cell_adt_qc_metrics", &per_cell_adt_qc_metrics);

    // suggest_adt_qc_filters.hpp
    mod.method("suggest_adt_qc_filters", &suggest_adt_qc_filters);

    // per_cell_crispr_qc_metrics.hpp
    mod.method("per_cell_crispr_qc_metrics", &per_cell_crispr_qc_metrics);

    // suggest_crispr_qc_filters.hpp
    mod.method("suggest_crispr_qc_filters", &suggest_crispr_qc_filters);

    // filter_cells.hpp
    mod.method("filter_cells", &filter_cells);
}
