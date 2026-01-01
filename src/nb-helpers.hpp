#pragma once
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/complex.h>
#include <vector>

namespace nb = nanobind;

// ZERO-COPY: Move a C++ vector into a NumPy array via a capsule
template <typename T>
auto MoveVectorToNdarray(std::vector<T>&& src) {
    auto* v_ptr = new std::vector<T>(std::move(src));
    size_t shape[1] = { v_ptr->size() };
    
    // The capsule acts as the owner, deleting the vector when Python GC hits
    nb::capsule owner(v_ptr, [](void *p) noexcept {
        delete reinterpret_cast<std::vector<T>*>(p);
    });

    return nb::ndarray<nb::numpy, T, nb::shape<-1>>(v_ptr->data(), 1, shape, owner);
}

// Optimization: Ensure input is contiguous and on CPU to avoid kernel logic errors
template <typename T>
std::vector<T> NdarrayToVector(nb::ndarray<T, nb::c_contig> arr) {
    return std::vector<T>(arr.data(), arr.data() + arr.size());
}

// Helper to move 2D vector to ndarray (flattened copy for now, or zero-copy if we change internal storage)
// Since internal storage is vector<vector>, we must copy to a flat buffer for numpy.
// To achieve zero-copy for 2D, we would need a flat vector in the C++ class.
// For now, we do a copy into a flat vector and then hand over ownership.
template <typename T>
auto MoveVector2DToNdarray(const std::vector<std::vector<T>>& src) {
    size_t rows = src.size();
    size_t cols = rows > 0 ? src[0].size() : 0;
    
    auto* flattened = new std::vector<T>();
    flattened->reserve(rows * cols);
    for (const auto& row : src) {
        flattened->insert(flattened->end(), row.begin(), row.end());
    }

    size_t shape[2] = { rows, cols };
    nb::capsule owner(flattened, [](void *p) noexcept { delete reinterpret_cast<std::vector<T>*>(p); });
    return nb::ndarray<nb::numpy, T, nb::shape<-1, -1>>(flattened->data(), 2, shape, owner);
}

// --- OUTPUT: Move nested vector to 2D Ndarray (Flattened Zero-Copy) ---
// Used for S1 and S2 batch results: [num_spheres][num_angles]
template <typename T>
auto MoveNestedVector2DToNdarray(std::vector<std::vector<T>>&& src) {
    size_t rows = src.size();
    size_t cols = (rows > 0) ? src[0].size() : 0;
    
    auto* flattened = new std::vector<T>();
    flattened->reserve(rows * cols);
    for (auto& row : src) {
        flattened->insert(flattened->end(), row.begin(), row.end());
    }

    size_t shape[2] = { rows, cols };
    nb::capsule owner(flattened, [](void *p) noexcept {
        delete reinterpret_cast<std::vector<T>*>(p);
    });

    return nb::ndarray<nb::numpy, T, nb::shape<-1, -1>>(
        flattened->data(), 2, shape, owner
    );
}

