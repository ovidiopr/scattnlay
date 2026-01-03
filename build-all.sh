#!/bin/bash
set -e

echo "--- 1. Installing Python Extension + CLI Tools ---"
# --no-build-isolation: Uses the libraries already in your environment (numpy, etc.)
# This prevents pip from creating a slow virtual environment for every build.
# STRICT_BUILD=ON ensures we don't skip SIMD/MP during dev.
python3 -m pip install -e . \
    -C build-dir=build \
    -C cmake.define.CMAKE_BUILD_TYPE=Release \
    -C cmake.define.WITH_HWY=ON \
    -C cmake.define.ENABLE_MP=ON
    
echo "--- 2. Building WASM Assets ---"
if [ -n "$EMSDK" ]; then
    echo "EMSDK found. Building WASM..."
    cmake --preset wasm-release -B build_wasm
    cmake --build build_wasm
else
    echo "WARNING: EMSDK not found. Skipping WASM build."
fi

echo "--- 3. Building Vue3 Frontend ---"
# WASM must still be handled separately as it uses a different toolchain (Emscripten)
cd guiapp
pnpm install
pnpm build
cd ..

echo "Build All completed."
