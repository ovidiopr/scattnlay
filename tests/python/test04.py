#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2017 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#
#    This file is part of scattnlay
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    The only additional remark is that we expect that all publications
#    describing work using this software, or all commercial products
#    using it, cite at least one of the following references:
#    [1] O. Peña and U. Pal, "Scattering of electromagnetic radiation by
#        a multilayered sphere," Computer Physics Communications,
#        vol. 180, Nov. 2009, pp. 2348-2354.
#    [2] K. Ladutenko, U. Pal, A. Rivera, and O. Peña-Rodríguez, "Mie
#        calculation of electromagnetic near-field for a multilayered
#        sphere," Computer Physics Communications, vol. 214, May 2017,
#        pp. 225-230.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

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

from scattnlay import scattnlay
import numpy as np

nL = 500.0
Xmax = 60.0

x = np.array([np.arange(1.0, nL + 1.0)*Xmax/nL], dtype = np.float64)
m = np.array([np.sqrt((2.0 - ((x[0] - 0.5*Xmax/nL)/60.0)**2.0)) + 0.0j], dtype = np.complex128)

theta = np.arange(0.0, 180.25, 0.25, dtype = np.float64)*np.pi/180.0

terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(x, m, theta)

S11 = S1[0].real*S1[0].real + S1[0].imag*S1[0].imag + S2[0].real*S2[0].real + S2[0].imag*S2[0].imag
result = np.vstack((theta*180.0/np.pi, S11/(2.0*Xmax*Xmax), np.cos(theta))).transpose()

try:
    import matplotlib.pyplot as plt

    plt.plot(result[ : , 0], result[ : , 1], 'k', result[ : , 0], result[ : , 2], 'r')

    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 1e3)

    plt.xlabel('Theta')

    plt.draw()
    plt.show()
finally:
    np.savetxt("test04.txt", result, fmt = "%.5f")
    print result


