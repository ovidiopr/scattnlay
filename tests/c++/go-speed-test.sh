#!/bin/bash
path=$PWD
PROGRAM='scattnlay.bin'

echo Compile with gcc
rm -f $PROGRAM

#google profiler  ######## FAST!!!
# file=speed-test.cc
# g++ -Ofast -std=c++11 $file ../../src/nmie.cc  -lm -lrt -o $PROGRAM /usr/lib/libtcmalloc.so.4 -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -march=native -mtune=native -msse4.2

file=speed-test.cc
#g++ -Ofast -std=c++11 $file ../../src/nmie.cc -DMULTI_PRECISION=20  -lm -lrt -o $PROGRAM /usr/lib/libtcmalloc_minimal.so.4 -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -march=native -mtune=native -msse4.2
g++ -Ofast -std=c++11 $file ../../src/nmie.cc  -lm -lrt -o $PROGRAM /usr/lib/libtcmalloc_minimal.so.4 -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -march=native -mtune=native -msse4.2

#file=speed-test-applied.cc
# g++ -Ofast -std=c++11 $file ../../src/nmie.cc ../../src/nmie-applied.cc -lm -lrt -o $PROGRAM /usr/lib/libtcmalloc_minimal.so.4 -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -march=native -mtune=native -msse4.2

echo Should be:
echo test01, +1.41154e+00, +4.17695e-01, +9.93844e-01, +1.59427e-01, +1.25809e+00, +3.67376e-01, +2.95915e-01
echo Result:
time  ./$PROGRAM -l 5 0.4642 1.8000 1.7000 0.7114 0.8000 0.7000 0.7393 1.2000 0.0900 0.9168 2.8000 0.2000 1.0000 1.5000 0.4000  -t 0.0 90.0 5 -c test01
# #result
