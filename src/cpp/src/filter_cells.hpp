#ifndef SCRANJL_FILTER_CELLS_HPP
#define SCRANJL_FILTER_CELLS_HPP

#include "jlcxx/jlcxx.hpp"
#include "ScranMatrix.h"
#include <cstdint>

#ifndef JL_BINDINGS_ONLY
#include "scran/scran.hpp"
#endif

ScranMatrix filter_cells(const ScranMatrix& x, jlcxx::ArrayRef<uint8_t> discard) 
#ifndef JL_BINDINGS_ONLY
{
    scran::FilterCells runner;
    return ScranMatrix(runner.run(x.ptr, discard.data()));
}
#endif
;

#endif
