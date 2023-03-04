#ifndef SCRANJL_MODEL_GENE_VAR_HPP
#define SCRANJL_MODEL_GENE_VAR_HPP

#include "jlcxx/jlcxx.hpp"
#include "ScranMatrix.h"
#include <cstdint>
#include <iostream>

#ifndef JL_BINDINGS_ONLY
#include "scran/scran.hpp"
#endif

struct ModelGeneVarResults {
    ModelGeneVarResults(scran::ModelGeneVar::Results res) : results(std::move(res)) {}

    jlcxx::ArrayRef<double> means(int b) {
        auto& m = results.means[b];
        return jlcxx::ArrayRef<double>(m.data(), m.size());
    }

    jlcxx::ArrayRef<double> variances(int b) {
        auto& m = results.variances[b];
        return jlcxx::ArrayRef<double>(m.data(), m.size());
    }

    jlcxx::ArrayRef<double> residuals(int b) {
        auto& m = results.residuals[b];
        return jlcxx::ArrayRef<double>(m.data(), m.size());
    }

    jlcxx::ArrayRef<double> fitted(int b) {
        auto& m = results.fitted[b];
        return jlcxx::ArrayRef<double>(m.data(), m.size());
    }

    scran::ModelGeneVar::Results results;
};

ModelGeneVarResults model_gene_var(const ScranMatrix& x, bool use_batch, jlcxx::ArrayRef<int32_t> batch, double span, int nthreads) 
#ifndef JL_BINDINGS_ONLY
{
    scran::ModelGeneVar mvar;
    mvar.set_num_threads(nthreads);
    mvar.set_span(span);
    const int32_t* bptr = (use_batch ? batch.data() : NULL);
    return ModelGeneVarResults(mvar.run_blocked(x.ptr.get(), bptr));
}
#endif
;

#endif
