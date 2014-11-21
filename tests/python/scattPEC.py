#!/usr/bin/env python

# This test case calculates the differential scattering
# cross section for different x values of a PEC sphere

# The differential cross section from wave optics is:
# d(Csca)/d(a**2*Omega) = S11(Theta)/x**2

from scattnlay import scattnlay
import numpy as np

dX = 0.5
Xmax = 5.0

m = np.array([[1.0 - 1.0j]], dtype = np.complex128)
theta = np.arange(0.0, 180.25, 0.25, dtype = np.float64)*np.pi/180.0

result = theta*180.0/np.pi

for xl in np.arange(dX, Xmax, dX, dtype = np.float64):
    x = np.array([[xl]], dtype = np.float64)
    terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(x, m, theta)

    S11 = S1[0].real*S1[0].real + S1[0].imag*S1[0].imag + S2[0].real*S2[0].real + S2[0].imag*S2[0].imag
    result = np.vstack((result, S11/(2.0*xl*xl)))

result = result.transpose()

try:
    import matplotlib.pyplot as plt

    plt.plot(result[ : , 0], result[ : , 1:])

    ax = plt.gca()
    ax.set_yscale('log')
#    ax.set_ylim(1e-4, 1e3)

    plt.xlabel('Theta')

    plt.draw()
    plt.show()
finally:
    np.savetxt("scattPEC.txt", result, fmt = "%.5f")
    print result


