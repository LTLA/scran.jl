#include "jlcxx/jlcxx.hpp"

#include "ScranMatrix.h"
#include "scran/scran.hpp"

#include <vector>
#include <cstdint>
#include <limits>
#include <type_traits>

void per_cell_rna_qc_metrics(const ScranMatrix& x, 
    jlcxx::ArrayRef<jl_value_t*> subsets, 
    jlcxx::ArrayRef<double> sums, 
    jlcxx::ArrayRef<int32_t> detected, 
    jlcxx::ArrayRef<jl_value_t*> proportions,
    int nthreads)
{
    std::vector<uint8_t*> subptrs;
    for (auto ptr : subsets) {
        jl_array_t* arr = reinterpret_cast<jl_array_t*>(ptr);
        subptrs.push_back(reinterpret_cast<uint8_t*>(jl_array_data(arr)));
    }

    scran::PerCellRnaQcMetrics::Buffers buffer;
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
