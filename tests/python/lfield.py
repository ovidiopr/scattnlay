#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2015 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#
#    This file is part of python-scattnlay
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
#    using it, cite the following reference:
#    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by
#        a multilayered sphere," Computer Physics Communications,
#        vol. 180, Nov. 2009, pp. 2348-2354.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# This test case calculates the electric field along the 
# X, Y and Z axes, for an spherical silver nanoparticle
# embedded in glass.

# Refractive index values correspond to a wavelength of
# 400 nm. Maximum of the surface plasmon resonance (and,
# hence, of electric field) is expected under those
# conditions.

from scattnlay import fieldnlay
import numpy as np

x = np.ones((1, 2), dtype = np.float64)
x[0, 0] = 2.0*np.pi*0.05/1.064
x[0, 1] = 2.0*np.pi*0.06/1.064

m = np.ones((1, 2), dtype = np.complex128)
m[0, 0] = 1.53413/1.3205
m[0, 1] = (0.565838 + 7.23262j)/1.3205

nc = 1001

coordX = np.zeros((nc, 3), dtype = np.float64)
coordY = np.zeros((nc, 3), dtype = np.float64)
coordZ = np.zeros((nc, 3), dtype = np.float64)

scan = np.linspace(-4.0*x[0, 1], 4.0*x[0, 1], nc)
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
    np.savetxt("lfield.txt", result, fmt = "%.5f")
    print result



