#!/bin/bash
echo Compile with gcc -O2
rm -rf *.bin
g++ -O2 standalone.cc nmie.cc ucomplex.cc -lm -o scattnlay.bin
cp scattnlay.bin ../scattnlay
cd tests/shell
# for file in `ls *.sh`; do ./$file; done
repeats=30
echo Run test for $repeats times
time for i in `seq $repeats`; do ./test01.sh; done
echo Run test with original binary for $repeats times
cp ../../../scattnlay-0.3.0 ../../../scattnlay
time for i in `seq $repeats`; do ./test01.sh; done
