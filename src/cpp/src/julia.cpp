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
#include "model_gene_var.hpp"
#include "run_pca.cpp"
#include "run_blocked_pca.cpp"

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
    mod.add_type<RunPcaResults>("RunPcaResults")
        .method("num_dims", &RunPcaResults::num_dims)
        .method("num_obs", &RunPcaResults::num_obs)
        .method("num_genes", &RunPcaResults::num_genes)
        .method("principal_components", &RunPcaResults::principal_components)
        .method("rotation", &RunPcaResults::rotation)
        .method("variance_explained", &RunPcaResults::variance_explained)
        .method("total_variance", &RunPcaResults::total_variance);

    mod.method("run_pca", &run_pca);

    // run_blocked_pca.hpp
    mod.add_type<BlockedPcaResults>("BlockedPcaResults")
        .method("num_dims", &BlockedPcaResults::num_dims)
        .method("num_obs", &BlockedPcaResults::num_obs)
        .method("num_genes", &BlockedPcaResults::num_genes)
        .method("principal_components", &BlockedPcaResults::principal_components)
        .method("rotation", &BlockedPcaResults::rotation)
        .method("variance_explained", &BlockedPcaResults::variance_explained)
        .method("total_variance", &BlockedPcaResults::total_variance);

    mod.method("run_blocked_pca", &run_blocked_pca);
}
