#include "suggest_rna_qc_filters.h"
#include <cstdint>

SuggestRnaQcFiltersResults suggest_rna_qc_filters(
    jlcxx::ArrayRef<double> sums, 
    jlcxx::ArrayRef<int32_t> detected, 
    jlcxx::ArrayRef<jl_value_t*> proportions, 
    bool use_block, 
    jlcxx::ArrayRef<int32_t> block_ids, 
    jlcxx::ArrayRef<uint8_t> keep, 
    double nmads) 
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

