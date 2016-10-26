#!/bin/bash
path=$PWD
file=compare.cc

echo Compile with gcc
rm -f *.bin
rm -f ../scattnlay

#g++ -Ofast -std=c++11 compare.cc nmie.cc  nmie-wrapper.cc -lm -lrt -o scattnlay.bin -static
# g++ -Ofast -std=c++11 compare.cc nmie.cc  nmie-wrapper.cc -lm -lrt -o scattnlay-pg.bin -static -pg

#google profiler  ######## FAST!!!
echo Uncomment next line to compile compare.cc
# g++ -Ofast -std=c++11 $file ../../src/nmie.cc  -lm -lrt -o $PROGRAM /usr/lib/libtcmalloc_minimal.so.4 -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -march=native -mtune=native -msse4.2

g++ -Ofast -std=c++11 $file nmie.cc  nmie-old.cc -lm -lrt -o scattnlay.bin /usr/lib/libtcmalloc.so.4 -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -march=native -mtune=native -msse4.2

#  g++ -Ofast -std=c++11 compare.cc nmie.cc  nmie-wrapper.cc -lm -lrt -o scattnlay-g.bin -ltcmalloc -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -g

#DEBUG!
#clang++ -g -O1 -fsanitize=address  -fno-optimize-sibling-calls -fno-omit-frame-pointer -std=c++11 compare.cc nmie.cc  nmie-old.cc -lm -lrt -o scattnlay.bin

cp scattnlay.bin ../scattnlay
# cp scattnlay-g.bin ../scattnlay-g
#cd tests/shell
# for file in `ls *.sh`;  do 
#     if [ "$file" != "test03.sh" ]; then
# 	./$file > /dev/null
# 	echo $file
#     fi
# done

PROGRAM='./scattnlay.bin'
time  ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.4  $PROGRAM -l 5 0.4642 1.8000 1.7000 0.7114 0.8000 0.7000 0.7393 1.2000 0.0900 0.9168 2.8000 0.2000 1.0000 1.5000 0.4000  -t 0.0 90.0 5 -c test01
# #result
# # test01, +1.41154e+00, +4.17695e-01, +9.93844e-01, +1.59427e-01, +1.25809e+00, +3.67376e-01, +2.95915e-01

# echo BUG -- All designs should give almost the same answer
# #echo
# echo $PROGRAM -l 1 4.71238898038469 2 0.0001
# $PROGRAM -l 1 4.71238898038469 2 0.0001
# #echo
# echo $PROGRAM -l 2 4.71238898038469 2 0.0001 9.42477796076937 1 0
# $PROGRAM -l 2 4.71238898038469 2 0.0001 9.42477796076937 1.5 0.0001
#echo
#echo $PROGRAM -l 2 4.71238898038469 2 0.0001 9.42477796076938 1 0
#$PROGRAM -l 2 4.71238898038469 2 0.0001 9.42477796076938 1.5 0.0001
#echo
#  #apt-get install oprofile
# echo oprofile
# PROGRAM='../../../scattnlay-g'
# rm -rf oprofiletmp
# rm -rf oprofile_data
# time operf $PROGRAM -l 5 0.4642 1.8000 1.7000 0.7114 0.8000 0.7000 0.7393 1.2000 0.0900 0.9168 2.8000 0.2000 1.0000 1.5000 0.4000  -t 0.0 90.0 5 -c test01
# #opreport --symbols > opreport.log
# mkdir oprofiletmp
# opannotate --source --output-dir=./oprofiletmp/

#echo valgring
# valgrind --tool=callgrind $PROGRAM -l 5 0.4642 1.8000 1.7000 0.7114 0.8000 0.7000 0.7393 1.2000 0.0900 0.9168 2.8000 0.2000 1.0000 1.5000 0.4000  -t 0.0 90.0 5 -c test01
# rm out.dot
# ./gprof2dot.py --output=out.dot --format=callgrind callgrind.out.*
# mv callgrind.out.* callgrind
# dot -Tsvg out.dot -o graph.svg

# rm *.aprof
# /home/mmedia/soft/aprof-0.2.1/inst/bin/valgrind --tool=aprof  $PROGRAM -l 5 0.4642 1.8000 1.7000 0.7114 0.8000 0.7000 0.7393 1.2000 0.0900 0.9168 2.8000 0.2000 1.0000 1.5000 0.4000  -t 0.0 90.0 5 -c test01


# rm gmon.out
# rm analysis.txt
# PROGRAM='../../../scattnlay-pg'
# time $PROGRAM -l 5 0.4642 1.8000 1.7000 0.7114 0.8000 0.7000 0.7393 1.2000 0.0900 0.9168 2.8000 0.2000 1.0000 1.5000 0.4000  -t 0.0 90.0 5 -c test01
# gprof $PROGRAM gmon.out > analysis.txt

# repeats=30
# echo Run test for $repeats times
# time for i in `seq $repeats`; do ./test01.sh; done
# echo Run test with original binary for $repeats times
# cp ../../../scattnlay-0.3.0 ../../../scattnlay
# time for i in `seq $repeats`; do ./test01.sh; done

#Run static analysis
# echo Run manually:
# echo make clean
# echo ~/../mmedia/soft/cov-analysis-linux64-7.6.0/bin/cov-build --dir cov-int make standalone
# echo      tar -zcvf cov-int.tar.gz cov-int
# echo to submit to coverity

# cd $path
# file=test-negative-epsilon.cc
# #g++ -Ofast -std=c++11 $file nmie.cc -lm -lrt -o $file.bin /usr/lib/libtcmalloc.so.4 -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -march=native -mtune=native -msse4.2

# #DEBUG!
# clang++ -g -O1 -fsanitize=address  -fno-optimize-sibling-calls -fno-omit-frame-pointer -std=c++11 $file nmie.cc -lm -lrt -o $file.bin

# ./$file.bin

### Cython
# rm scattnlay.so
# export CFLAGS='-std=c++11'
# python setup.py build_ext --inplace
# cp scattnlay.so tests/python/
# cd tests/python/
#./field-Ag-flow.py
# ./lfield.py
# ./field-dielectric-sphere.py
# ./field.py
# ./test01.py
# ./test01-wrapper.py
