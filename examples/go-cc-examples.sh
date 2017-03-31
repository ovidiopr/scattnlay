#!/bin/bash
path=$PWD
PROGRAM='scattnlay-example.bin'

file=example-eval-force.cc
echo Compile $file with gcc
rm -f $PROGRAM
g++ -Ofast -std=c++11 $file  ../src/shell-generator.cc -lm -lrt -o $PROGRAM -march=native -mtune=native -msse4.2 -Wall
# g++ -Ofast -std=c++11 $file ../src/nmie.cc ../src/nmie-applied.cc ../src/shell-generator.cc -lm -lrt -o $PROGRAM -march=native -mtune=native -msse4.2
echo Compilation done. Running...
time ./$PROGRAM

# file=test-surf-integral.cc
# echo Compile $file with gcc
# rm -f $PROGRAM
# g++ -Ofast -std=c++11 $file  ../src/shell-generator.cc -lm -lrt -o $PROGRAM -march=native -mtune=native -msse4.2
# # g++ -Ofast -std=c++11 $file ../src/nmie.cc ../src/nmie-applied.cc ../src/shell-generator.cc -lm -lrt -o $PROGRAM -march=native -mtune=native -msse4.2
# echo Compilation done. Running...
# ./$PROGRAM


# file=example-minimal.cc
# echo Compile $file with gcc
# rm -f $PROGRAM
# g++ -Ofast -std=c++11 $file ../src/nmie.cc ../src/nmie-applied.cc -lm -lrt -o $PROGRAM -march=native -mtune=native -msse4.2
# echo Compilation done. Running...
# ./$PROGRAM


# file=example-get-Mie.cc
# echo Compile $file with gcc
# rm -f $PROGRAM
# # Production for multiprecision
# #g++ -Ofast -std=c++11 $file ../src/nmie.cc ../src/nmie-applied.cc read-spectra.cc -DMULTI_PRECISION=200 -lm -lrt -o $PROGRAM /usr/lib/libtcmalloc_minimal.so.4 -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -march=native -mtune=native -msse4.2

# # Simplified for multiprecision
# #g++ -Ofast -std=c++11 $file ../src/nmie.cc ../src/nmie-applied.cc ./read-spectra.cc -DMULTI_PRECISION=200 -lm -lrt -o $PROGRAM -march=native -mtune=native -msse4.2

# # Simplified for double precision
# g++ -Ofast -std=c++11 $file ../src/nmie.cc ../src/nmie-applied.cc ./read-spectra.cc -lm -lrt -o $PROGRAM -march=native -mtune=native -msse4.2
# echo Compilation done. Running...
# ./$PROGRAM
