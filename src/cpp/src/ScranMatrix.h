#ifndef SCRAN_MATRIX_H
#define SCRAN_MATRIX_H

#include "tatami/tatami.hpp"
#include "jlcxx/jlcxx.hpp"

struct ScranMatrix {
    ScranMatrix(tatami::NumericMatrix* p) : ptr(p) {}

    ScranMatrix(std::shared_ptr<tatami::NumericMatrix> p) : ptr(std::move(p)) {}

    int nrow() const {
        return ptr->nrow();
    }

    int ncol() const {
        return ptr->ncol();
    }

    void row(size_t r, jlcxx::ArrayRef<double> buffer) const {
        ptr->row_copy(r, buffer.data());
        return;
    }

    void column(size_t c, jlcxx::ArrayRef<double> buffer) const {
        ptr->column_copy(c, buffer.data());
        return;
    }

    std::shared_ptr<tatami::NumericMatrix> ptr;
};

#endif
