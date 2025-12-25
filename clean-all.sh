#!/bin/bash

echo "Cleaning build directories..."
rm -rf build build_std build_mp build_simd build_wasm build_deps .tox _deps

echo "Cleaning Python artifacts..."
rm -rf build/ dist/ # Python setup.py build dir and dist
rm -rf scattnlay.egg-info
find . -name "*.so" -delete
find . -name "*.dylib" -delete
find . -name "*.pyd" -delete
find . -name "__pycache__" -type d -exec rm -rf {} +

echo "Cleaning guiapp artifacts..."
rm -f guiapp/src/nmiejs.js
rm -f guiapp/src/nmiejs.wasm

echo "Cleaning other artifacts..."
rm -f scattnlay.bin scattnlay-pg.bin # From go.sh if used
rm -f *.bin

echo "Clean complete."
