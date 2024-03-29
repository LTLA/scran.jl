#include "jlcxx/jlcxx.hpp"

#define JL_BINDINGS_ONLY
#include "initialize_from_memory.hpp"
#include "log_norm_counts.hpp"
#include "grouped_size_factors.hpp"
#include "per_cell_rna_qc_metrics.hpp"
#include "suggest_rna_qc_filters.hpp"
#include "per_cell_adt_qc_metrics.hpp"
#include "suggest_adt_qc_filters.hpp"
#include "per_cell_crispr_qc_metrics.hpp"
#include "suggest_crispr_qc_filters.hpp"
#include "filter_cells.hpp"
#include "model_gene_var.hpp"
#include "run_pca.cpp"
#include "run_blocked_pca.cpp"
#include "run_multibatch_pca.cpp"

template<class Results>
void register_pca_results_class(jlcxx::Module& mod, const std::string& name) {
    mod.add_type<Results>(name)
        .method("num_dims", &Results::num_dims)
        .method("num_obs", &Results::num_obs)
        .method("num_genes", &Results::num_genes)
        .method("principal_components", &Results::principal_components)
        .method("rotation", &Results::rotation)
        .method("variance_explained", &Results::variance_explained)
        .method("total_variance", &Results::total_variance);
}

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

    mod.add_type<ModelGeneVarResults>("ModelGeneVarResults")
        .method("means", &ModelGeneVarResults::means)
        .method("variances", &ModelGeneVarResults::variances)
        .method("residuals", &ModelGeneVarResults::residuals)
        .method("fitted", &ModelGeneVarResults::fitted);

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

    // grouped_size_factors.hpp
    mod.method("grouped_size_factors", &grouped_size_factors);

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

    // model_gene_var.hpp
    mod.method("model_gene_var", &model_gene_var);

    // run_pca.hpp
    register_pca_results_class<RunPcaResults>(mod, "RunPcaResults");

    mod.method("run_pca", &run_pca);

    // run_blocked_pca.hpp
    register_pca_results_class<BlockedPcaResults>(mod, "BlockedPcaResults");

    mod.method("run_blocked_pca", &run_blocked_pca);

    // run_multibatch_pca.hpp
    register_pca_results_class<MultiBatchPcaResults>(mod, "MultiBatchPcaResults");

    mod.method("run_multibatch_pca", &run_multibatch_pca);
}
