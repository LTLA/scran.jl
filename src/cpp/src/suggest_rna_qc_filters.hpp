#ifndef SUGGEST_RNA_QC_FILTERS_H
#define SUGGEST_RNA_QC_FILTERS_H

#include "jlcxx/jlcxx.hpp"
#include "scran/scran.hpp"
#include <cstdint>

struct SuggestRnaQcFiltersResults {
    SuggestRnaQcFiltersResults(scran::SuggestRnaQcFilters::Thresholds r) : results(std::move(r)) {}

    jlcxx::ArrayRef<double> thresholds_sums() {
        return jlcxx::ArrayRef<double>(results.sums.data(), results.sums.size());
    }

    jlcxx::ArrayRef<double> thresholds_detected() {
        return jlcxx::ArrayRef<double>(results.detected.data(), results.detected.size());
    }

    jlcxx::ArrayRef<double> thresholds_proportions(int i) {
        return jlcxx::ArrayRef<double>(results.subset_proportions[i].data(), results.subset_proportions[i].size());
    }

    scran::SuggestRnaQcFilters::Thresholds results;
};

SuggestRnaQcFiltersResults suggest_rna_qc_filters(
    jlcxx::ArrayRef<double> sums, 
    jlcxx::ArrayRef<int32_t> detected, 
    jlcxx::ArrayRef<jl_value_t*> proportions, 
    bool use_block, 
    jlcxx::ArrayRef<int32_t> block_ids, 
    jlcxx::ArrayRef<uint8_t> keep, 
    double nmads) 
#ifndef JL_BINDINGS_ONLY
{
    scran::PerCellRnaQcMetrics::Buffers<double, int32_t> buffer;
    buffer.sums = sums.data();
    buffer.detected = detected.data();
    for (auto ptr : proportions) {
        jl_array_t* arr = reinterpret_cast<jl_array_t*>(ptr);
        buffer.subset_proportions.push_back(reinterpret_cast<double*>(jl_array_data(arr)));
    }

    scran::SuggestRnaQcFilters runner;
    runner.set_num_mads(nmads);

    const int32_t* bptr = (use_block ? block_ids.data() : NULL);
    auto res = runner.run_blocked(sums.size(), bptr, buffer);
    res.filter_blocked(sums.size(), bptr, buffer, keep.data());

    return SuggestRnaQcFiltersResults(std::move(res));
}
#endif
;

#endif
