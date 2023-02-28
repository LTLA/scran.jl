#ifndef SUGGEST_RNA_QC_FILTERS_H
#define SUGGEST_RNA_QC_FILTERS_H

#include "jlcxx/jlcxx.hpp"
#include "scran/scran.hpp"

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

#endif
