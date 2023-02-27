#include "jlcxx/jlcxx.hpp"

#include "ScranMatrix.h"

#include <vector>
#include <cstdint>
#include <limits>
#include <type_traits>

template<class XVector, class IVector, class PVector>
ScranMatrix create_matrix_copy_byrow(XVector x, IVector i, PVector p, int nrow, int ncol, bool byrow) {
    for (auto& idx : i) {
        --idx; // get back to 0-based.
    }

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

template<typename DataType, class RowType>
ScranMatrix create_matrix_no_copy(const jlcxx::ArrayRef<DataType>& x, std::vector<RowType> i, std::vector<size_t> p, int nrow, int ncol, bool byrow) {
    tatami::ArrayView<DataType> x_(x.data(), x.size());
    for (auto& idx : i) {
        --idx; // get back to 0-based.
    }

    if (byrow) {
        typedef tatami::CompressedSparseMatrix<true, double, int, decltype(x_), decltype(i), decltype(p)> SparseMat;
        return ScranMatrix(new SparseMat(nrow, ncol, std::move(x_), std::move(i), std::move(p), false));
    } else {
        typedef tatami::CompressedSparseMatrix<false, double, int, decltype(x_), decltype(i), decltype(p)> SparseMat;
        return ScranMatrix(new SparseMat(nrow, ncol, std::move(x_), std::move(i), std::move(p), false));
    }
}

template<typename DataType, typename IndexType, typename PtrType>
ScranMatrix initialize_from_memory(jlcxx::ArrayRef<DataType> x, jlcxx::ArrayRef<IndexType> i, jlcxx::ArrayRef<PtrType> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced) {
    std::vector<size_t> p_(p.begin(), p.end());
    for (auto& pidx : p_) {
        --pidx; // get back to 0-based.
    }

    // If no copy is requested, we hold onto the array views.
    if (no_copy) {
        if (nrow <= std::numeric_limits<uint16_t>::max()) {
            return create_matrix_no_copy(x, std::vector<uint16_t>(i.begin(), i.end()), std::move(p_), nrow, ncol, byrow);
        } else {
            return create_matrix_no_copy(x, std::vector<int>(i.begin(), i.end()), std::move(p_), nrow, ncol, byrow);
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

ScranMatrix initialize_from_memory_int32_int32(jlcxx::ArrayRef<int32_t> x, jlcxx::ArrayRef<int32_t> i, jlcxx::ArrayRef<int32_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced) {
    return initialize_from_memory(std::move(x), std::move(i), std::move(p), nrow, ncol, no_copy, byrow, forced);
}

ScranMatrix initialize_from_memory_int64_int32(jlcxx::ArrayRef<int64_t> x, jlcxx::ArrayRef<int32_t> i, jlcxx::ArrayRef<int32_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced) {
    return initialize_from_memory(std::move(x), std::move(i), std::move(p), nrow, ncol, no_copy, byrow, forced);
}

ScranMatrix initialize_from_memory_float32_int32(jlcxx::ArrayRef<float> x, jlcxx::ArrayRef<int32_t> i, jlcxx::ArrayRef<int32_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced) {
    return initialize_from_memory(std::move(x), std::move(i), std::move(p), nrow, ncol, no_copy, byrow, forced);
}

ScranMatrix initialize_from_memory_float64_int32(jlcxx::ArrayRef<double> x, jlcxx::ArrayRef<int32_t> i, jlcxx::ArrayRef<int32_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced) {
    return initialize_from_memory(std::move(x), std::move(i), std::move(p), nrow, ncol, no_copy, byrow, forced);
}

ScranMatrix initialize_from_memory_int32_int64(jlcxx::ArrayRef<int32_t> x, jlcxx::ArrayRef<int64_t> i, jlcxx::ArrayRef<int64_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced) {
    return initialize_from_memory(std::move(x), std::move(i), std::move(p), nrow, ncol, no_copy, byrow, forced);
}

ScranMatrix initialize_from_memory_int64_int64(jlcxx::ArrayRef<int32_t> x, jlcxx::ArrayRef<int64_t> i, jlcxx::ArrayRef<int64_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced) {
    return initialize_from_memory(std::move(x), std::move(i), std::move(p), nrow, ncol, no_copy, byrow, forced);
}

ScranMatrix initialize_from_memory_float32_int64(jlcxx::ArrayRef<float> x, jlcxx::ArrayRef<int64_t> i, jlcxx::ArrayRef<int64_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced) {
    return initialize_from_memory(std::move(x), std::move(i), std::move(p), nrow, ncol, no_copy, byrow, forced);
}

ScranMatrix initialize_from_memory_float64_int64(jlcxx::ArrayRef<double> x, jlcxx::ArrayRef<int64_t> i, jlcxx::ArrayRef<int64_t> p, int nrow, int ncol, bool no_copy, bool byrow, bool forced) {
    return initialize_from_memory(std::move(x), std::move(i), std::move(p), nrow, ncol, no_copy, byrow, forced);
}
