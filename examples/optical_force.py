#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2017 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#    Copyright (C) 2013-2017 Konstantin Ladutenko <kostyfisik@gmail.com>
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

# This test case calculates the optical force over a silver nanoparticle,
# as a function of the irradiance and the radius.

from scattnlay import scattnlay
import numpy as np
from scipy.constants import pi, c

radius = np.linspace(0.5, 180.0, 360)
nAg = np.sqrt(-4.0 + 0.7j)
wl = 400.0

x = np.ones((len(radius), 1), dtype = np.float64)
x[:, 0] = 2.0*pi*radius/wl

m = np.ones((len(radius), 1), dtype = np.complex128)
m[:, 0] *= nAg

terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(x, m)
F = pi*Qpr*radius*radius/c/1e9

result = np.vstack((radius, 1e11*F, 1e13*F, 1e15*F)).transpose()

try:
    import matplotlib.pyplot as plt

    plt.figure(1)
    plt.subplot(311)
    plt.plot(radius, 1e11*F, 'k', label = '10$^{11}$ W/m$^2$')
    plt.plot(radius, 1e13*F, 'b', label = '10$^{13}$ W/m$^2$')
    plt.plot(radius, 1e15*F, 'g', label = '10$^{15}$ W/m$^2$')
    plt.ylabel('F (nN)')
    plt.legend(loc = 4)
    ax = plt.gca()
    ax.set_yscale('log')

    plt.subplot(312)
    plt.plot(radius, g, 'r', label = 'g')
    plt.ylabel('g')

    plt.subplot(313)
    plt.plot(radius, Qext, 'k', label = 'Q$_{ext}$')
    plt.plot(radius, Qsca, 'b', label = 'Q$_{sca}$')
    plt.plot(radius, Qpr, 'g', label = 'Q$_{pr}$')
    plt.ylabel('Q')
    plt.legend()

    plt.xlabel('R (nm)')
    
    plt.show()
finally:
    #np.savetxt("test_force.txt", result, fmt = "%.5e")
    print result

