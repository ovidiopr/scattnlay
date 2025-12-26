#!/bin/bash
set -e


echo "Testing..."
if [ -d "build" ]; then
    cd build
    ctest --output-on-failure
    cd ..
else
    echo "build not found, skipping."
fi

echo "Running Vitest..."
if [ -d "guiapp" ]; then
    cd guiapp
    pnpm install
    pnpm test
    cd ..
fi

echo "All tests completed."
