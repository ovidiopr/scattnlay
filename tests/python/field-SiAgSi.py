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
# E-k plane, for an spherical Si-Ag-Si nanoparticle. Core radius is 17.74 nm,
# inner layer 23.31nm, outer layer 22.95nm. Working wavelength is 800nm, we use
# silicon epsilon=13.64+i0.047, silver epsilon= -28.05+i1.525

import scattnlay
from scattnlay import fieldnlay
from scattnlay import scattnlay
import numpy as np

epsilon_Si = 13.64 + 0.047j
epsilon_Ag = -28.05 + 1.525j
# epsilon_Si = 2.0 + 0.047j
# epsilon_Ag = -2.0 + 1.525j

index_Si = np.sqrt(epsilon_Si)
index_Ag = np.sqrt(epsilon_Ag)

# Values for 800 nm, taken from http://refractiveindex.info/
index_Si = 3.69410 + 0.0065435j
index_Ag = 0.18599 + 4.9886j

WL=800 #nm
core_width = 17.74 #nm Si
inner_width = 23.31 #nm Ag
outer_width = 22.95 #nm  Si

core_r = core_width
inner_r = core_r+inner_width
outer_r = inner_r+outer_width

# n1 = 1.53413
# n2 = 0.565838 + 7.23262j
nm = 1.0

x = np.ones((1, 3), dtype = np.float64)
x[0, 0] = 2.0*np.pi*core_r/WL
x[0, 1] = 2.0*np.pi*inner_r/WL
x[0, 2] = 2.0*np.pi*outer_r/WL

m = np.ones((1, 3), dtype = np.complex128)
m[0, 0] = index_Si/nm
m[0, 1] = index_Ag/nm
m[0, 2] = index_Si/nm

print "x =", x
print "m =", m

npts = 281

scan = np.linspace(-2.0*x[0, 2], 2.0*x[0, 2], npts)

coordX, coordZ = np.meshgrid(scan, scan)
coordX.resize(npts*npts)
coordZ.resize(npts*npts)
coordY = np.zeros(npts*npts, dtype = np.float64)

coord = np.vstack((coordX, coordY, coordZ)).transpose()

terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(x, m)
terms, E, H = fieldnlay(x, m, coord)
print("Qabs = "+str(Qabs));
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

    edata = np.resize(Eh, (npts, npts)).T

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Rescale to better show the axes
    scale_x = np.linspace(min(coordX)*1.064/2.0/np.pi/nm, max(coordX)*1.064/2.0/np.pi/nm, npts)
    scale_z = np.linspace(min(coordZ)*1.064/2.0/np.pi/nm, max(coordZ)*1.064/2.0/np.pi/nm, npts)

    # Define scale ticks
    min_tick = min(min_tick, np.amin(edata))
    max_tick = max(max_tick, np.amax(edata))
    # scale_ticks = np.power(10.0, np.linspace(np.log10(min_tick), np.log10(max_tick), 6))
    scale_ticks = np.linspace(min_tick, max_tick, 6)

    # Interpolation can be 'nearest', 'bilinear' or 'bicubic'
    cax = ax.imshow(edata, interpolation = 'nearest', cmap = cm.jet,
                    origin = 'lower', vmin = min_tick, vmax = max_tick,
                    extent = (min(scale_x), max(scale_x), min(scale_z), max(scale_z))
                    #,norm = LogNorm()
                    )

    # Add colorbar
    cbar = fig.colorbar(cax, ticks = [a for a in scale_ticks])
    cbar.ax.set_yticklabels(['%3.1e' % (a) for a in scale_ticks]) # vertically oriented colorbar
    pos = list(cbar.ax.get_position().bounds)
    fig.text(pos[0] - 0.02, 0.925, '|E|/|E$_0$|', fontsize = 14)

    plt.xlabel('Z')
    plt.ylabel('X')

    # # This part draws the nanoshell
    # from matplotlib import patches

    # s1 = patches.Arc((0, 0), 1.0*x[0, 0], 1.0*x[0, 0],  angle=0.0, zorder=2,
    #                  theta1=0.0, theta2=360.0, linewidth=1, color='black')
    # ax.add_patch(s1)

    # s2 = patches.Arc((0, 0), 2.0*x[0, 1], 2.0*x[0, 1], angle=0.0, zorder=2,
    #                  theta1=0.0, theta2=360.0, linewidth=1, color='black')
    # ax.add_patch(s2) 
    # s3 = patches.Arc((0, 0), 2.0*x[0, 2], 2.0*x[0, 2], angle=0.0, zorder=2,
    #                  theta1=0.0, theta2=360.0, linewidth=1, color='black')
    # ax.add_patch(s3) 
    # # End of drawing

    plt.draw()

    plt.show()

    plt.clf()
    plt.close()
finally:
    np.savetxt("field.txt", result, fmt = "%.5f")
    print result


