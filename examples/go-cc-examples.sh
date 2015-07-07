#!/bin/bash
path=$PWD
PROGRAM='scattnlay-example.bin'

echo Compile with gcc
rm -f $PROGRAM

file=example-get-Mie.cc
g++ -Ofast -std=c++11 $file ../src/nmie.cc ../src/nmie-applied.cc -lm -lrt -o $PROGRAM /usr/lib/libtcmalloc.so.4 -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -march=native -mtune=native -msse4.2

./$PROGRAM
# #result
