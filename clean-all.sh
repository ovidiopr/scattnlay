#!/bin/bash
# 1. Clean C++ Core (Shared Build Directory)
rm -rf build/ build_wasm/ build_native/

# 2. Clean Python artifacts
rm -rf .tox/ dist/ *.egg-info/
find scattnlay -name "*.so" -delete
find scattnlay -name "*.pyd" -delete
find . -name "__pycache__" -type d -exec rm -rf {} +

# 3. Clean Vue3 GUI
cd guiapp
rm -rf dist/ 
rm public/wasm/nmiejs.*
cd ..

echo "Environment cleaned."
