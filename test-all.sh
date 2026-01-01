#!/bin/bash
set -e

echo "--- Running C++ Tests (Includes Tox via CTest) ---"
# Since your CMakeLists.txt adds tox as a test, this covers C++ and Python.
ctest --test-dir build --output-on-failure

echo "--- Running Vue3/WASM Tests ---"
cd guiapp
pnpm test
cd ..

echo "All tests passed."
