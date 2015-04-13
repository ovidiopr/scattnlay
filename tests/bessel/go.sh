#!/bin/bash
clang++ -g -O1 -fsanitize=address  -fno-optimize-sibling-calls -fno-omit-frame-pointer -std=c++11 ../../bessel.cc test_bessel.cc -lm -lrt -o test_bessel.bin
./test_bessel.bin
