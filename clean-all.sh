#!/bin/bash
echo "--- Cleaning C++ Core ---"
rm -rf build/ build_debug/ build_simd/ build_mp/ build_wasm/

echo "--- Cleaning Vue3 Frontend ---"
rm -rf guiapp/dist/ guiapp/node_modules/

echo "--- Cleaning Python Extension ---"
rm -rf .tox/ dist/ scattnlay.egg-info/
find scattnlay -name "*.so" -delete
find scattnlay -name "*.pyd" -delete

echo "Cleanup complete."
