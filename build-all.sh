#!/bin/bash
set -e

echo "Building SIMD (Google Highway) and Multi-Precision..."
cmake -S . -B build -DWITH_HWY=ON -DENABLE_MP=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

if command -v emcmake >/dev/null 2>&1; then
    echo "Building WebAssembly..."
    emcmake cmake -S . -B build_wasm -DCMAKE_BUILD_TYPE=Release
    cmake --build build_wasm -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
else
    echo "emcmake not found, skipping WebAssembly build."
fi

echo "All builds completed successfully."
