rm -rf build_native && python3 -m pip install -e . \
    -C build-dir=build_native \
    -C cmake.define.CMAKE_BUILD_TYPE=Release \
    -C cmake.define.WITH_HWY=ON \
    -C cmake.define.ENABLE_MP=ON
python3 examples/calc-simd-benchmark-nearfield.py
python3 examples/calc-simd-benchmark.py