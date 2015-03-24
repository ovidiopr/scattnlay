#!/usr/bin/env python

# This test case calculates the differential scattering
# cross section from a Luneburg lens, as described in:
# B. R. Johnson, Applied Optics 35 (1996) 3286-3296.

# The Luneburg lens is a sphere of radius a, with a
# radially-varying index of refraction, given by:
# m(r) = [2 - (r/a)**1]**(1/2)

# For the calculations, the Luneburg lens was approximated
# as a multilayered sphere with 500 equally spaced layers.
# The refractive index of each layer is defined to be equal to
# m(r) at the midpoint of the layer: ml = [2 - (xm/xL)**1]**(1/2),
# with xm = (xl-1 + xl)/2, for l = 1,2,...,L. The size
# parameter in the lth layer is xl = l*xL/500. According to
# geometrical optics theory, the differential cross section
# can be expressed as:
# d(Csca)/d(a**2*Omega) = cos(Theta)

# The differential cross section from wave optics is:
# d(Csca)/d(a**2*Omega) = S11(Theta)/x**2

from scattnlay import fieldnlay
import numpy as np

x = np.ones((1, 1), dtype = np.float64)
x[0, 0] = 1.

m = np.ones((1, 1), dtype = np.complex128)
m[0, 0] = (0.0252 + 2.0181j)/1.46

nc = 1001

coordX = np.zeros((nc, 3), dtype = np.float64)
coordY = np.zeros((nc, 3), dtype = np.float64)
coordZ = np.zeros((nc, 3), dtype = np.float64)

scan = np.linspace(-10.0*x[0, 0], 10.0*x[0, 0], nc)
one = np.ones(nc, dtype = np.float64)

coordX[:, 0] = scan
coordY[:, 1] = scan
coordZ[:, 2] = scan

terms, Ex, Hx = fieldnlay(x, m, coordX)
terms, Ey, Hy = fieldnlay(x, m, coordY)
terms, Ez, Hz = fieldnlay(x, m, coordZ)

Exr = np.absolute(Ex)
Eyr = np.absolute(Ey)
Ezr = np.absolute(Ez)

# |E|/|Eo|
Exh = np.sqrt(Exr[0, :, 0]**2 + Exr[0, :, 1]**2 + Exr[0, :, 2]**2)
Eyh = np.sqrt(Eyr[0, :, 0]**2 + Eyr[0, :, 1]**2 + Eyr[0, :, 2]**2)
Ezh = np.sqrt(Ezr[0, :, 0]**2 + Ezr[0, :, 1]**2 + Ezr[0, :, 2]**2)

result = np.vstack((scan, Exh, Eyh, Ezh)).transpose()

try:
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.errorbar(result[:, 0], one, fmt = 'k')
    ax.errorbar(result[:, 0], result[:, 1], fmt = 'r', label = 'X axis')
    ax.errorbar(result[:, 0], result[:, 2], fmt = 'g', label = 'Y axis')
    ax.errorbar(result[:, 0], result[:, 3], fmt = 'b', label = 'Z axis')

    ax.legend()

    plt.xlabel('X|Y|Z')
    plt.ylabel('|E|/|Eo|')

    plt.draw()
    plt.show()

finally:
    np.savetxt("field.txt", result, fmt = "%.5f")
    print result


