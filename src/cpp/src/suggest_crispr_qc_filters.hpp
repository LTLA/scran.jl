#ifndef SCRANJL_SUGGEST_CRISPR_QC_FILTERS_H
#define SCRANJL_SUGGEST_CRISPR_QC_FILTERS_H

#include "jlcxx/jlcxx.hpp"
#include "scran/scran.hpp"
#include <cstdint>

struct SuggestCrisprQcFiltersResults {
    SuggestCrisprQcFiltersResults(scran::SuggestCrisprQcFilters::Thresholds r) : results(std::move(r)) {}

    jlcxx::ArrayRef<double> thresholds_max_count() {
        return jlcxx::ArrayRef<double>(results.max_count.data(), results.max_count.size());
    }

    scran::SuggestCrisprQcFilters::Thresholds results;
};

SuggestCrisprQcFiltersResults suggest_crispr_qc_filters(
    jlcxx::ArrayRef<double> sums, 
    jlcxx::ArrayRef<double> max_proportion, 
    bool use_block, 
    jlcxx::ArrayRef<int32_t> block_ids, 
    jlcxx::ArrayRef<uint8_t> keep, 
    double nmads) 
#ifndef JL_BINDINGS_ONLY
{
    scran::PerCellCrisprQcMetrics::Buffers<double, int32_t> buffer;
    buffer.sums = sums.data();
    buffer.max_proportion = max_proportion.data();

    scran::SuggestCrisprQcFilters runner;
    runner.set_num_mads(nmads);

    const int32_t* bptr = (use_block ? block_ids.data() : NULL);
    auto res = runner.run_blocked(sums.size(), bptr, buffer);
    res.filter_blocked(sums.size(), bptr, buffer, keep.data());

    return SuggestCrisprQcFiltersResults(std::move(res));
}
#endif
;

#endif
