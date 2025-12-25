#!/bin/bash
set -e


echo "Building Multi-Precision..."
cmake -S . -B build_mp -DENABLE_MP=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_mp -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

echo "Building SIMD (Google Highway)..."
cmake -S . -B build_simd -DWITH_HWY=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_simd -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

if command -v emcmake >/dev/null 2>&1; then
    echo "Building WebAssembly..."
    emcmake cmake -S . -B build_wasm -DCMAKE_BUILD_TYPE=Release
    cmake --build build_wasm -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
else
    echo "emcmake not found, skipping WebAssembly build."
fi

echo "All builds completed successfully."
