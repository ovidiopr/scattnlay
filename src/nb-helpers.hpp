#pragma once
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/complex.h>
#include <vector>
#include <complex>

namespace nb = nanobind;

// 1D: Move C++ vector to NumPy without copy
template <typename T>
auto MoveVectorToNdarray(std::vector<T>&& src) {
    auto* v_ptr = new std::vector<T>(std::move(src));
    size_t shape[1] = { v_ptr->size() };
    nb::capsule owner(v_ptr, [](void *p) noexcept { delete reinterpret_cast<std::vector<T>*>(p); });
    return nb::ndarray<nb::numpy, T, nb::shape<-1>>(v_ptr->data(), 1, shape, owner);
}

// 2D: Flatten vector<vector> into a contiguous block for NumPy
template <typename T>
auto MoveVector2DToNdarray(std::vector<std::vector<T>>&& src) {
    size_t rows = src.size();
    size_t cols = (rows > 0) ? src[0].size() : 0;
    auto* flattened = new std::vector<T>();
    flattened->reserve(rows * cols);
    for (auto& row : src) {
        flattened->insert(flattened->end(), std::make_move_iterator(row.begin()), std::make_move_iterator(row.end()));
    }
    size_t shape[2] = { rows, cols };
    nb::capsule owner(flattened, [](void *p) noexcept { delete reinterpret_cast<std::vector<T>*>(p); });
    return nb::ndarray<nb::numpy, T, nb::shape<-1, -1>>(flattened->data(), 2, shape, owner);
}

// Optimization: Ensure input is contiguous and on CPU to avoid kernel logic errors
template <typename T>
std::vector<T> NdarrayToVector(nb::ndarray<T, nb::c_contig> arr) {
    return std::vector<T>(arr.data(), arr.data() + arr.size());
}

// Conversion helpers
template <typename TargetT, typename SourceT>
std::vector<TargetT> ConvertVector(const std::vector<SourceT>& src) {
    std::vector<TargetT> dst;
    dst.reserve(src.size());
    for (const auto& val : src) {
        dst.push_back(static_cast<TargetT>(val));
    }
    return dst;
}

template <typename TargetT, typename SourceT>
std::vector<std::complex<TargetT>> ConvertComplexVector(const std::vector<std::complex<SourceT>>& src) {
    std::vector<std::complex<TargetT>> dst;
    dst.reserve(src.size());
    for (const auto& val : src) {
        dst.emplace_back(static_cast<TargetT>(val.real()), static_cast<TargetT>(val.imag()));
    }
    return dst;
}

template <typename TargetT, typename SourceT>
std::vector<std::vector<std::complex<TargetT>>> ConvertComplexVectorVector(const std::vector<std::vector<std::complex<SourceT>>>& src) {
    std::vector<std::vector<std::complex<TargetT>>> dst;
    dst.reserve(src.size());
    for (const auto& row : src) {
        dst.push_back(ConvertComplexVector<TargetT>(row));
    }
    return dst;
}

