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

echo "All tests completed."
