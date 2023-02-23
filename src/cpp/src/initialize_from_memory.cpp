#include "jlcxx/jlcxx.hpp"

#include "tatami/tatami.hpp"

#include <vector>
#include <cstdint>
#include <limits>
#include <type_traits>

struct ScranMatrix {
    ScranMatrix(tatami::NumericMatrix* p) : ptr(p) {}

    void row(size_t r, jlcxx::ArrayRef<double> buffer) {
        ptr->row_copy(r, buffer.data());
        return;
    }

    std::shared_ptr<tatami::NumericMatrix> ptr;
};

template<class XVector, class IVector, class PVector>
ScranMatrix create_matrix_copy_byrow(XVector x, IVector i, PVector p, int nrow, int ncol, bool byrow) {
    if (byrow) {
        typedef tatami::CompressedSparseMatrix<true, double, int, XVector, IVector, PVector> SparseMat;
        return ScranMatrix(new SparseMat(nrow, ncol, std::move(x), std::move(i), std::move(p), false));
    } else {
        typedef tatami::CompressedSparseMatrix<false, double, int, XVector, IVector, PVector> SparseMat;
        return ScranMatrix(new SparseMat(nrow, ncol, std::move(x), std::move(i), std::move(p), false));
    }
}

template<class Incoming, class RowType>
ScranMatrix create_matrix_copy_x_max(const Incoming& x, std::vector<RowType> i, std::vector<size_t> p, int nrow, int ncol, bool byrow) {
    auto mined = (x.size() ? *std::min_element(x.begin(), x.end()) : 0);
    if (mined < 0) {
        throw std::runtime_error("expression values should be positive");
    }

    auto maxed = (x.size() ? *std::max_element(x.begin(), x.end()) : 0);
    if (maxed <= std::numeric_limits<uint16_t>::max()) {
        return create_matrix_copy_byrow(std::vector<uint16_t>(x.begin(), x.end()), std::move(i), std::move(p), nrow, ncol, byrow);
    } else {
        return create_matrix_copy_byrow(std::vector<int>(x.begin(), x.end()), std::move(i), std::move(p), nrow, ncol, byrow);
    }
}

template<typename DataType, class RowType>
ScranMatrix create_matrix_copy_x_type(const jlcxx::ArrayRef<DataType>& x, std::vector<RowType> i, std::vector<size_t> p, int nrow, int ncol, bool byrow, bool forced) {
    if constexpr(std::is_same<DataType, int>::value) {
        return create_matrix_copy_x_max(x, std::move(i), std::move(p), nrow, ncol, byrow);
    } else {
        if (forced) {
            return create_matrix_copy_x_max(x, std::move(i), std::move(p), nrow, ncol, byrow);
        } else {
            return create_matrix_copy_byrow(std::vector<double>(x.begin(), x.end()), std::move(i), std::move(p), nrow, ncol, byrow);
        }
    }
}

template<typename DataType>
ScranMatrix initialize_from_memory(jlcxx::ArrayRef<DataType> x, jlcxx::ArrayRef<int> i, jlcxx::ArrayRef<int64_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced) {
    std::vector<size_t> p_(p.begin(), p.end());

    // If no copy is requested, we hold onto the array views.
    if (no_copy) {
        tatami::ArrayView<DataType> x_(x.data(), x.size());
        tatami::ArrayView<int> i_(i.data(), i.size());

        if (byrow) {
            typedef tatami::CompressedSparseMatrix<true, double, int, decltype(x_), decltype(i_), decltype(p)> SparseMat;
            return ScranMatrix(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p), false));
        } else {
            typedef tatami::CompressedSparseMatrix<false, double, int, decltype(x_), decltype(i_), decltype(p)> SparseMat;
            return ScranMatrix(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p), false));
        }
    }

    // Otherwise, beginning the copying process with the indices.
    if (nrow <= std::numeric_limits<uint16_t>::max()) {
        return create_matrix_copy_x_type(x, std::vector<uint16_t>(i.begin(), i.end()), std::move(p_), nrow, ncol, byrow, forced);
    } else {
        return create_matrix_copy_x_type(x, std::vector<int>(i.begin(), i.end()), std::move(p_), nrow, ncol, byrow, forced);
    }
}

/*** Defining concrete bindings. ***/

ScranMatrix initialize_from_memory_integer(jlcxx::ArrayRef<int> x, jlcxx::ArrayRef<int> i, jlcxx::ArrayRef<int64_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced) {
    return initialize_from_memory(std::move(x), std::move(i), std::move(p), nrow, ncol, no_copy, byrow, forced);
}

ScranMatrix initialize_from_memory_double(jlcxx::ArrayRef<double> x, jlcxx::ArrayRef<int> i, jlcxx::ArrayRef<int64_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced) {
    return initialize_from_memory(std::move(x), std::move(i), std::move(p), nrow, ncol, no_copy, byrow, forced);
}

JLCXX_MODULE define_module_initialize_from_memory(jlcxx::Module& mod) {
    mod.add_type<ScranMatrix>("ScranMatrix")
       .method("row", &ScranMatrix::row);

    mod.method("initialize_from_memory_integer", &initialize_from_memory_integer);
    mod.method("initialize_from_memory_double", &initialize_from_memory_double);
}


