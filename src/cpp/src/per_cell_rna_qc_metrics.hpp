#ifndef SCRANJL_PER_CELL_RNA_QC_METRICS_HPP
#define SCRANJL_PER_CELL_RNA_QC_METRICS_HPP

#include "jlcxx/jlcxx.hpp"
#include "ScranMatrix.h"
#include <cstdint>

#ifndef JL_BINDINGS_ONLY
#include "scran/scran.hpp"
#include <vector>
#endif

void per_cell_rna_qc_metrics(const ScranMatrix& x, 
    jlcxx::ArrayRef<jl_value_t*> subsets, 
    jlcxx::ArrayRef<double> sums, 
    jlcxx::ArrayRef<int32_t> detected, 
    jlcxx::ArrayRef<jl_value_t*> proportions,
    int nthreads)
#ifndef JL_BINDINGS_ONLY
{
    std::vector<uint8_t*> subptrs;
    for (auto ptr : subsets) {
        jl_array_t* arr = reinterpret_cast<jl_array_t*>(ptr);
        subptrs.push_back(reinterpret_cast<uint8_t*>(jl_array_data(arr)));
    }

    scran::PerCellRnaQcMetrics::Buffers<double, int32_t> buffer;
    buffer.sums = sums.data();
    buffer.detected = detected.data();
    for (auto ptr : proportions) {
        jl_array_t* arr = reinterpret_cast<jl_array_t*>(ptr);
        buffer.subset_proportions.push_back(reinterpret_cast<double*>(jl_array_data(arr)));
    }

    scran::PerCellRnaQcMetrics runner;
    runner.set_num_threads(nthreads);
    runner.run(x.ptr.get(), subptrs, buffer);
    return;
}
#endif
;

#endif
