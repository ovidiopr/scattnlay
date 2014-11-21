#!/bin/bash
# This test case calculates the differential scattering
# cross section from a Luneburg lens, as described in:
# B. R. Johnson, Applied Optics 35 (1996) 3286-3296.
#
# The Luneburg lens is a sphere of radius a, with a
# radially-varying index of refraction, given by:
# m(r) = [2 - (r/a)**1]**(1/2)
#
# For the calculations, the Luneburg lens was approximated
# as a multilayered sphere with 500 equally spaced layers.
# The refractive index of each layer is defined to be equal to
# m(r) at the midpoint of the layer: ml = [2 - (xm/xL)**1]**(1/2),
# with xm = (xl-1 + xl)/2, for l = 1,2,...,L. The size
# parameter in the lth layer is xl = l*xL/500. According to
# geometrical optics theory, the differential cross section
# can be expressed as:
# d(Csca)/d(a**2*Omega) = cos(Theta)
#
# The differential cross section from wave optics is:
# d(Csca)/d(a**2*Omega) = S11(Theta)/x**2

PROGRAM='../../../scattnlay -l'

L=500
xL=60

FULL_PROGRAM="$PROGRAM $L"
for l in `seq 1 $L`; # Total number of layers
do
  xlM1=$(echo "scale=2; ($l-1)*$xL/$L" | bc -l)
  xl=$(echo "scale=2; $l*$xL/$L" | bc -l)
  xm=$(echo "scale=2; ($xlM1+$xl)/2" | bc -l)
  ml=$(echo "scale=5; e(0.5*l(2-($xm/$xL)*($xm/$xL)))" | bc -l)

  FULL_PROGRAM="$FULL_PROGRAM $xl $ml 0.0"
done

FULL_PROGRAM="$FULL_PROGRAM -t 0.0 180.0 721"
$FULL_PROGRAM

