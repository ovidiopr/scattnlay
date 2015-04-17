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

# This test case calculates the electric field in the 
# XY plane, for an spherical silver nanoparticle
# embedded in glass.

# Refractive index values correspond to a wavelength of
# 400 nm. Maximum of the surface plasmon resonance (and,
# hence, of electric field) is expected under those
# conditions.
import scattnlay

import os
path = os.path.dirname(scattnlay.__file__)
print(scattnlay.__file__)

from scattnlay import fieldnlay
import numpy as np

x = np.ones((1, 1), dtype = np.float64)
x[0, 0] = 0.001
#x[0, 1] = 2.0*np.pi*0.06/1.064

m = np.ones((1, 1), dtype = np.complex128)
m[0, 0] = 4.0
#m[0, 1] = (0.565838 + 7.23262j)/1.3205

npts = 501

scan = np.linspace(-2.0*x[0, 0], 2.0*x[0, 0], npts)

coordX, coordY = np.meshgrid(scan, scan)
coordX.resize(npts*npts)
coordY.resize(npts*npts)
coordZ = np.zeros(npts*npts, dtype = np.float64)

coord = np.vstack((coordX, coordY, coordZ)).transpose()

terms, E, H = fieldnlay(x, m, coord)

Er = np.absolute(E)

# |E|/|Eo|
#Eh = np.sqrt(Er[0, :, 0]**2 + Er[0, :, 1]**2 + Er[0, :, 2]**2)
Eh = Er[0, :, 0]**2 + Er[0, :, 1]**2 + Er[0, :, 2]**2

result = np.vstack((coordX, coordY, coordZ, Eh)).transpose()

try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import LogNorm

    min_tick = 0.25
    max_tick = 1.0

    edata = np.resize(Eh, (npts, npts))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Rescale to better show the axes
    scale_x = np.linspace(min(coordX), max(coordX), npts)
    scale_y = np.linspace(min(coordY), max(coordY), npts)

    # Define scale ticks
    min_tick = max(0.1, min(min_tick, np.amin(edata)))
    max_tick = max(max_tick, np.amax(edata))
    #scale_ticks = np.power(10.0, np.linspace(np.log10(min_tick), np.log10(max_tick), 6))
    scale_ticks = np.linspace(min_tick,max_tick, 5)
    #scale_ticks = np.linspace(0, 2, 11)

    # Interpolation can be 'nearest', 'bilinear' or 'bicubic'
    cax = ax.imshow(edata, interpolation = 'nearest', cmap = cm.jet,
                    origin = 'lower', vmin = min_tick, vmax = max_tick,
                    #origin = 'lower', vmin = 0.16, vmax = 0.18,
                    extent = (min(scale_x), max(scale_x), min(scale_y), max(scale_y))
                    #,norm = LogNorm()
                    )

    # Add colorbar
    cbar = fig.colorbar(cax, ticks = [a for a in scale_ticks])
    cbar.ax.set_yticklabels(['%3.1e' % (a) for a in scale_ticks]) # vertically oriented colorbar
    pos = list(cbar.ax.get_position().bounds)
    fig.text(pos[0] - 0.02, 0.925, '|E|/|E$_0$|', fontsize = 14)

    plt.xlabel('X')
    plt.ylabel('Y')

    # # This part draws the nanoshell
    # from matplotlib import patches

    # s1 = patches.Arc((0, 0), 2.0*x[0, 0], 2.0*x[0, 0], angle=0.0, zorder=2,
    #                   theta1=0.0, theta2=360.0, linewidth=1, color='#00fa9a')
    # ax.add_patch(s1)

    # s2 = patches.Arc((0, 0), 2.0*x[0, 0], 2.0*x[0, 0], angle=0.0, zorder=2,
    #                   theta1=0.0, theta2=360.0, linewidth=1, color='#00fa9a')
    # ax.add_patch(s2)
    # # End of drawing

    plt.draw()

    plt.show()

    plt.clf()
    plt.close()
finally:
    np.savetxt("field.txt", result, fmt = "%.5f")
    print result


