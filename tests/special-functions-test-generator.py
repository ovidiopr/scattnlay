#!/usr/bin/env python3
import mpmath as mp
import numpy as np
Du_test = [
# // x, [Re(m), Im(m)], Qext, Qsca, test_name
[0.099, [0.75,0], 7.417859e-06, 7.417859e-06, 'a'],
# [0.101, [0.75,0], 8.033538e-06, 8.033538e-06, 'b'],
# [10,    [0.75,0],     2.232265, 2.232265, 'c'],
# [1000,  [0.75,0],     1.997908, 1.997908, 'd'],
# [100,   [1.33,-1e-5], 2.101321, 2.096594, 'e'],
# [10000, [1.33,-1e-5], 2.004089, 1.723857, 'f'],
# [0.055, [1.5, -1],    0.101491, 1.131687e-05, 'g'],
# [0.056, [1.5, -1],   0.1033467, 1.216311e-05, 'h'],
# [100,   [1.5, -1],    2.097502, 1.283697, 'i'],
# [10000, [1.5, -1],    2.004368, 1.236574, 'j'],
# [1,     [10,  -10],   2.532993, 2.049405, 'k'],
# [100,   [10,  -10,],  2.071124, 1.836785, 'l'],
[10000, [10,  -10],   2.005914, 1.795393, 'm'],
[80, [1.05,  1],   0, 0, 'Yang'],
]
# // Dtest refractive index is m={1.05,1}, the size parameter is x = 80
n_list = [0,1,30,50,60,70,75,80,85,90,99,116,130];

def get_z_values(du_list):
    zlist = []
    for record in du_list:
        x = mp.mpf(str(record[0]))
        m = mp.mpc(str(record[1][0]), str(record[1][1]))
        z = x*m
        zlist.append(z)
    return zlist

# Riccati-Bessel z*j_n(z)
def psi(n,z):
    return mp.sqrt( (mp.pi * z)/2 ) * mp.besselj(n+1/2,z)

# to compare r(n,z) with Wolfram Alpha
# n=49, z=1.3-2.1i,  SphericalBesselJ[n-1,z]/SphericalBesselJ[n,z]
def r(n,z):
    return psi(n-1,z)/psi(n,z)

def D1(n,z):
    return r(n,z) - n/z

def main ():
    mp.mp.dps = 62
    output_dps = 10

    zlist = get_z_values(Du_test)
    for z in zlist:
        print(z)
        for n in n_list:
            print('{',mp.nstr(D1(n,z).real, output_dps),',',
                  mp.nstr(D1(n,z).imag, output_dps),'}')

main()
