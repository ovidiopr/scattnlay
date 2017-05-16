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

# This test case calculates the the electric field in the
# XY plane, for a Luneburg lens, as described in:
# B. R. Johnson, Applied Optics 35 (1996) 3286-3296.

# The Luneburg lens is a sphere of radius a, with a
# radially-varying index of refraction, given by:
# m(r) = [2 - (r/a)**1]**(1/2)

# For the calculations, the Luneburg lens was approximated
# as a multilayered sphere with 500 equally spaced layers.
# The refractive index of each layer is defined to be equal to
# m(r) at the midpoint of the layer: ml = [2 - (xm/xL)**1]**(1/2),
# with xm = (xl-1 + xl)/2, for l = 1,2,...,L. The size
# parameter in the lth layer is xl = l*xL/500.

from scattnlay import fieldnlay
import numpy as np

nL = 500.0
Xmax = 60.0

x = np.ones((1, nL), dtype = np.float64)
x[0] = np.arange(1.0, nL + 1.0)*Xmax/nL

m = np.ones((1, nL), dtype = np.complex128)
m[0] = np.sqrt((2.0 - ((x[0] - 0.5*Xmax/nL)/60.0)**2.0)) + 0.0j

print "x =", x
print "m =", m

npts = 501

scan = np.linspace(-10.0*x[0, -1], 10.0*x[0, -1], npts)

coordX, coordY = np.meshgrid(scan, scan)
coordX.resize(npts*npts)
coordY.resize(npts*npts)
coordZ = np.zeros(npts*npts, dtype = np.float64)

coord = np.vstack((coordX, coordY, coordZ)).transpose()

terms, E, H = fieldnlay(x, m, coord)

Er = np.absolute(E)

# |E|/|Eo|
Eh = np.sqrt(Er[0, :, 0]**2 + Er[0, :, 1]**2 + Er[0, :, 2]**2)

result = np.vstack((coordX, coordY, coordZ, Eh)).transpose()

try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import LogNorm

    min_tick = 0.1
    max_tick = 1.0

    edata = np.resize(Eh, (npts, npts))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Rescale to better show the axes
    scale_x = np.linspace(min(coordX), max(coordX), npts)
    scale_y = np.linspace(min(coordY), max(coordY), npts)

    # Define scale ticks
    min_tick = min(min_tick, np.amin(edata))
    max_tick = max(max_tick, np.amax(edata))
    scale_ticks = np.power(10.0, np.linspace(np.log10(min_tick), np.log10(max_tick), 6))

    # Interpolation can be 'nearest', 'bilinear' or 'bicubic'
    cax = ax.imshow(edata, interpolation = 'nearest', cmap = cm.jet,
                    origin = 'lower', vmin = min_tick, vmax = max_tick,
                    extent = (min(scale_x), max(scale_x), min(scale_y), max(scale_y)),
                    norm = LogNorm())

    # Add colorbar
    cbar = fig.colorbar(cax, ticks = [a for a in scale_ticks])
    cbar.ax.set_yticklabels(['%3.1e' % (a) for a in scale_ticks]) # vertically oriented colorbar
    pos = list(cbar.ax.get_position().bounds)
    fig.text(pos[0] - 0.02, 0.925, '|E|/|E$_0$|', fontsize = 14)

    plt.xlabel('X')
    plt.ylabel('Y')

    # This part draws the nanoshell
#    from matplotlib import patches

#    s1 = patches.Arc((0, 0), 2.0*x[0, 0], 2.0*x[0, 0], angle=0.0, zorder=2,
#                      theta1=0.0, theta2=360.0, linewidth=1, color='#00fa9a')
#    ax.add_patch(s1)

#    s2 = patches.Arc((0, 0), 2.0*x[0, 1], 2.0*x[0, 1], angle=0.0, zorder=2,
#                      theta1=0.0, theta2=360.0, linewidth=1, color='#00fa9a')
#    ax.add_patch(s2)
    # End of drawing

    plt.draw()

    plt.show()

    plt.clf()
    plt.close()
finally:
    np.savetxt("test04_field.txt", result, fmt = "%.5f")
    print result


