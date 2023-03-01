#ifndef SCRANJL_SUGGEST_ADT_QC_FILTERS_H
#define SCRANJL_SUGGEST_ADT_QC_FILTERS_H

#include "jlcxx/jlcxx.hpp"
#include "scran/scran.hpp"
#include <cstdint>

struct SuggestAdtQcFiltersResults {
    SuggestAdtQcFiltersResults(scran::SuggestAdtQcFilters::Thresholds r) : results(std::move(r)) {}

    jlcxx::ArrayRef<double> thresholds_detected() {
        return jlcxx::ArrayRef<double>(results.detected.data(), results.detected.size());
    }

    jlcxx::ArrayRef<double> thresholds_totals(int i) {
        return jlcxx::ArrayRef<double>(results.subset_totals[i].data(), results.subset_totals[i].size());
    }

    scran::SuggestAdtQcFilters::Thresholds results;
};

SuggestAdtQcFiltersResults suggest_adt_qc_filters(
    jlcxx::ArrayRef<int32_t> detected, 
    jlcxx::ArrayRef<jl_value_t*> totals, 
    bool use_block, 
    jlcxx::ArrayRef<int32_t> block_ids, 
    jlcxx::ArrayRef<uint8_t> keep, 
    double min_detected_drop,
    double nmads) 
#ifndef JL_BINDINGS_ONLY
{
    scran::PerCellAdtQcMetrics::Buffers<double, int32_t> buffer;
    buffer.detected = detected.data();
    for (auto ptr : totals) {
        jl_array_t* arr = reinterpret_cast<jl_array_t*>(ptr);
        buffer.subset_totals.push_back(reinterpret_cast<double*>(jl_array_data(arr)));
    }

    scran::SuggestAdtQcFilters runner;
    runner
        .set_num_mads(nmads)
        .set_min_detected_drop(min_detected_drop);

    const int32_t* bptr = (use_block ? block_ids.data() : NULL);
    auto res = runner.run_blocked(detected.size(), bptr, buffer);
    res.filter_blocked(detected.size(), bptr, buffer, keep.data());

    return SuggestAdtQcFiltersResults(std::move(res));
}
#endif
;

#endif
