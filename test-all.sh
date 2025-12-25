#!/bin/bash
set -e

echo "Testing Standard (Double Precision)..."
if [ -d "build_std" ]; then
    cd build_std
    ctest --output-on-failure
    cd ..
else
    echo "build_std not found, skipping."
fi

echo "Testing Multi-Precision..."
if [ -d "build_mp" ]; then
    cd build_mp
    ctest --output-on-failure
    cd ..
else
    echo "build_mp not found, skipping."
fi

echo "Testing SIMD..."
if [ -d "build_simd" ]; then
    cd build_simd
    ctest --output-on-failure
    cd ..
else
    echo "build_simd not found, skipping."
fi

echo "All tests completed."
