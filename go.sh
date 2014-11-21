#!/bin/bash
echo Run test of complex libs
cd ./tests/dev
rm -rf *.bin
clang++ test-complex-lib.cc ../../ucomplex.cc -o test-complex-lib.cc.bin
./test-complex-lib.cc.bin
