#!/usr/bin/env python3
import mpmath as mp


# APPLIED OPTICS / Vol. 53, No. 31 / 1 November 2014, eq(13)
def LeRu_cutoff(z):
    x = mp.fabs(z)
    return int(x + 11 * x**(1/3) + 1)


# Riccati-Bessel z*j_n(z)
def psi(n,z):
    return mp.sqrt( (mp.pi * z)/2 ) * mp.autoprec(mp.besselj)(n+1/2,z)


# to compare r(n,z) with Wolfram Alpha
# n=49, z=1.3-2.1i,  SphericalBesselJ[n-1,z]/SphericalBesselJ[n,z]
def r(n,z):
    if n > 0:
        return psi(n-1,z)/psi(n,z)
    return mp.cos(z)/mp.sin(z)


def D1(n,z):
    return r(n,z) - n/z
