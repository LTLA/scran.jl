#ifndef SCRANJL_PER_CELL_CRISPR_QC_METRICS_HPP
#define SCRANJL_PER_CELL_CRISPR_QC_METRICS_HPP

#include "jlcxx/jlcxx.hpp"
#include "ScranMatrix.h"
#include <cstdint>

#ifndef JL_BINDINGS_ONLY
#include "scran/scran.hpp"
#include <vector>
#endif

void per_cell_crispr_qc_metrics(const ScranMatrix& x, 
    jlcxx::ArrayRef<double> sums, 
    jlcxx::ArrayRef<int32_t> detected, 
    jlcxx::ArrayRef<double> max_proportion,
    jlcxx::ArrayRef<int32_t> max_index,
    int nthreads)
#ifndef JL_BINDINGS_ONLY
{
    scran::PerCellCrisprQcMetrics::Buffers<double, int32_t> buffer;
    buffer.sums = sums.data();
    buffer.detected = detected.data();
    buffer.max_proportion = max_proportion.data();
    buffer.max_index = max_index.data();

    scran::PerCellCrisprQcMetrics runner;
    runner.set_num_threads(nthreads);
    runner.run(x.ptr.get(), buffer);
    return;
}
#endif
;

#endif
