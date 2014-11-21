#!/bin/bash
# This test explores the applicability limits of scattnlay
# for Kai & Massoli's model:
# L. Kai and P. Massoli, Applied Optics 33 (1994) 501-511.
#
# The model consists in a multilayered sphere with a radial 
# profile of the refractive index ml = nl + i*kl given by:
# nl = n1 + 0.5(nL - n1)(1 - cos(t*Theta)) and kl = 0.0,
# where t = (l - 1)/(L - 1), n1 = 1.01nL, and nL = 1.33.
# The size parameter is xl = x1 + t(xL - x1), where
# l = 1,2,...,L, x1 = 0.001xL and L is the total number of
# layers.
#
# WARNING: This test can take a LOT of time to completely run.

PROGRAM='../../../scattnlay -l'

numL=200
maxL=1000
num_xL=800
max_xL=4000
nL=1.33

n1=$(echo "scale=5; 1.01*$nL" | bc -l)
pi=$(echo "scale=10; 4*a(1)" | bc -l)

for i in `seq 1 $numL`; # Total number of layers
do
  L=$(echo "scale=0; $i*$maxL/($numL)" | bc -l)
  for j in `seq 1 $num_xL`; # Total number of size factors
  do
    xL=$(echo "scale=2; $j*$max_xL/($num_xL)" | bc -l)
    x1=$(echo "scale=5; 0.001*$xL" | bc -l)
    FULL_PROGRAM="$PROGRAM $L"
    for l in `seq 1 $L`;
    do
      t=$(echo "scale=5; ($l-1)/($L-1)" | bc -l)
      n=$(echo "scale=2; $n1+0.5*($nL-$n1)*(1-c($t*$pi))" | bc -l)
      xl=$(echo "scale=2; $x1+$t*($xL-$x1)" | bc -l)
      FULL_PROGRAM="$FULL_PROGRAM $xl $n 0.0"
    done
    FULL_PROGRAM="$FULL_PROGRAM -c L=$L,xL=$xL"
    $FULL_PROGRAM
  done
done 
