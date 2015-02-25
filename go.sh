#!/bin/bash
echo Compile with gcc -O2
rm -rf *.bin

#g++ -O2 -std=c++11 standalone.cc nmie.cc -lm -o scattnlay.bin

clang++ -g -O1 -fsanitize=address  -fno-optimize-sibling-calls -fno-omit-frame-pointer -std=c++11 standalone.cc nmie.cc nmie-wrapper.cc -lm -o scattnlay.bin

cp scattnlay.bin ../scattnlay
cd tests/shell
# for file in `ls *.sh`;  do 
#     if [ "$file" != "test03.sh" ]; then
# 	./$file > /dev/null
# 	echo $file
#     fi
# done
PROGRAM='../../../scattnlay'
ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.4 $PROGRAM -l 5 0.4642 1.8000 1.7000 0.7114 0.8000 0.7000 0.7393 1.2000 0.0900 0.9168 2.8000 0.2000 1.0000 1.5000 0.4000  -t 0.0 90.0 5 -c test01


# repeats=30
# echo Run test for $repeats times
# time for i in `seq $repeats`; do ./test01.sh; done
# echo Run test with original binary for $repeats times
# cp ../../../scattnlay-0.3.0 ../../../scattnlay
# time for i in `seq $repeats`; do ./test01.sh; done
