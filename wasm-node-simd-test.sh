#!/bin/bash
set -e

# 1. Setup paths
PROJECT_ROOT=$(pwd)
BUILD_DIR="$PROJECT_ROOT/build_wasm"
DEPS_DIR="$PROJECT_ROOT/build_wasm_deps"
WASM_OUT_DIR="$PROJECT_ROOT/guiapp/public/wasm"

# 2. Check for Emscripten Environment
if [ -z "$EMSDK" ]; then
    echo "Error: EMSDK environment variable not found."
    echo "Please source your emsdk_env.sh before running this script."
    exit 1
fi

echo "--- Initializing Directories ---"
mkdir -p "$BUILD_DIR"
mkdir -p "$DEPS_DIR"
mkdir -p "$WASM_OUT_DIR"

# 3. Configure and Build
# We set FETCHCONTENT_BASE_DIR to our persistent deps folder.
# We use the toolchain file provided by the Emscripten SDK.
echo "--- Running CMake Configuration ---"
cmake -S . -B "$BUILD_DIR" \
    -DCMAKE_TOOLCHAIN_FILE="$EMSDK/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake" \
    -DCMAKE_BUILD_TYPE=Release \
    -DFETCHCONTENT_BASE_DIR="$DEPS_DIR" \
    -DWITH_HWY=ON \
    -DBUILD_EXAMPLES=OFF \
    -DBUILD_PYTHON_EXT=OFF

echo "--- Building nmiejs ---"
# Building specific target 'nmiejs'
cmake --build "$BUILD_DIR" --target nmiejs -j $(nproc || sysctl -n hw.ncpu)

# 4. Verify output exists
if [ ! -f "$WASM_OUT_DIR/nmiejs.js" ]; then
    echo "Error: Build failed to produce nmiejs.js in $WASM_OUT_DIR"
    exit 1
fi

# 5. Run Node.js SIMD Benchmark
echo "--- Running Node.js SIMD Benchmark ---"
# We go to the test directory to ensure relative imports work if needed, 
# or run directly pointing to the build artifacts.
cd "$PROJECT_ROOT/tests/wasm_node"

# Ensure node modules are installed for the test if package.json exists
if [ -f "package.json" ]; then
    npm install
fi

# Run the SIMD benchmark script using Node.js
# Note: Node needs --experimental-wasm-simd on older versions, 
# but modern Node (20+) supports it by default.
node wasm_simd_bench.js

echo "--- Test Complete ---"