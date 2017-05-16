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

# This is a test against the program n-mie (version 3a) for a water sphere surrounded by soot
# n-mie is based in the algorithm described in:
# Wu Z.P., Wang Y.P.
# Electromagnetic scattering for multilayered spheres:
# recursive algorithms
# Radio Science 1991. V. 26. P. 1393-1401.
# Voshchinnikov N.V., Mathis J.S.
# Calculating Cross Sections of Composite Interstellar Grains
# Astrophys. J. 1999. V. 526. #1. 

# The refractive indices of water and soot are m1 1.33 i0.00, m2 1.59 i0.66, respectively.
# The volume fraction of soot is 0.01.
# This test case was described in:
# W. Yang, Appl. Opt. 42 (2003) 1710-1720.

from scattnlay import scattnlay
import numpy as np

x = np.ones((400, 2), dtype = np.float64)
x[:, 1] = np.arange(0.1, 100.0, 0.25)
x[:, 0] = 0.99**(1.0/3.0)*x[:, 1]

m = np.ones((400, 2), dtype = np.complex128)
m[:, 0] *= 1.33
m[:, 1] *= 1.59 + 0.66j

terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(x, m)

result = np.vstack((x[:, 1], Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo)).transpose()

try:
    import matplotlib.pyplot as plt

    plt.figure(1)
    plt.subplot(311)
    plt.plot(x[:, 1], Qext, 'k')
    plt.ylabel('Qext')

    plt.subplot(312)
    plt.plot(x[:, 1], Qsca, 'r')
    plt.ylabel('Qsca')

    plt.subplot(313)
    plt.plot(x[:, 1], Albedo, 'g')
    plt.ylabel('Albedo')

    plt.xlabel('X')

    plt.show()
finally:
    np.savetxt("test02.txt", result, fmt = "%.5f")
    print result

