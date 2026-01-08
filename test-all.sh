#!/bin/bash
set -e

echo "--- 1. Running C++ Tests & Python Tox (via CTest) ---"
# This runs your GoogleTest binaries AND triggers Tox (as defined in your CMakeLists.txt)
ctest --test-dir build --output-on-failure

echo "--- 2. Running Vue3/Frontend Vitest ---"
cd guiapp
pnpm test
cd ..

tox

./build/tests/test_farfield_simd_benchmark | grep -e "speedup:" -e "time:" -e "avg" -e "RUN" -e "SIMD"
python3 tests/test_simd_benchmarks.py 
./build/tests/test_nearfield_simd_benchmark | grep -e "Speedup:" -e "time:" -e "Avg" -e "RUN"
echo "All tests passed!"