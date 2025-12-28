#!/bin/bash
clear
echo \`\`\`bash
NAME=scattnlay
echo RUN repo2txt
touch ${NAME}.txt
rm ${NAME}.txt && repo2txt -r . -o ${NAME}.txt --exclude-dir 'dist' 'refractiveindex.info-database' 'tests/shell' build build_wasm build_deps build_mp build_simd python_scattnlay.egg-info 'scattnlay/.tox' .tox .git __pycache__ build doc debian vue-cli3-webapp examples node_modules CMakeFiles build_bench Temporary --ignore-types .wasm .js .pyc --ignore-files scattnlay.txt test_near_field_multi_precision test_Riccati_Bessel_logarithmic_derivative gprof2dot.py LICENSE README.md opreport.log COPYING PKG-INFO test01.txt test02.txt test03.txt test04.txt test04_fin.txt nmiejs.wasm nmiejs.js fieldnlay-dp scattnlay-dp test_spec_functions_data.hpp nearfield-dp farfield-dp bench_special_functions test_bulk_sphere test_bulk_sphere_multi_precision test_near_field_multi_precision test_Riccati_Bessel_logarithmic_derivative  CMakeCache.txt Makefile
echo \`\`\`
