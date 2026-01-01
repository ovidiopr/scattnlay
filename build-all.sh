#!/bin/bash
set -e

echo "--- Building C++ Core (Default Preset) ---"
cmake --preset default
cmake --build --preset default

echo "--- Building Vue3 Frontend (Pnpm/Quasar) ---"
cd guiapp
pnpm install
pnpm build
cd ..

echo "--- Building Python Extension (Editable/Strict) ---"
# We use editable mode for development. STRICT_BUILD=ON ensures
# we have Boost and Highway properly configured.
pip install -e . -v -Ccmake.args="-DSTRICT_BUILD=ON"
