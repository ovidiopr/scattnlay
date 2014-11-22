#!/bin/bash
echo Run test of complex libs
echo **Debug
echo ***Clang
rm -rf *.bin
clang++ -std=c++11 test-complex-lib.cc ../../ucomplex.cc -o test-complex-lib.cc.bin
./test-complex-lib.cc.bin
echo ***Gcc
rm -rf *.bin
g++ -std=c++11 test-complex-lib.cc ../../ucomplex.cc -o test-complex-lib.cc.bin
./test-complex-lib.cc.bin
echo **Build
echo ***Clang
rm -rf *.bin
clang++ -std=c++11 -O2 test-complex-lib.cc ../../ucomplex.cc -o test-complex-lib.cc.bin
./test-complex-lib.cc.bin
echo ***Gcc
rm -rf *.bin
g++ -std=c++11 -O2  test-complex-lib.cc ../../ucomplex.cc -o test-complex-lib.cc.bin
./test-complex-lib.cc.bin
