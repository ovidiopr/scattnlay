#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2016 Paul Müller (paul.mueller [at] biotec.tu-dresden.de)
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
#    using it, cite the following reference:
#    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by
#        a multilayered sphere," Computer Physics Communications,
#        vol. 180, Nov. 2009, pp. 2348-2354.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# This test case calculates the phase retardation that is introduced
# by a weak dielectric sphere to an incident plane wave. Only the
# x-polarized light is considered.

# Note: This example computes the phase behind the sphere. In microscopy,
# the focal plane during imaging is close to the center of the sphere.
# To compute the phase image that corresponds to that imaged with a focal
# plane at the center of the bead, numerical refocusing of the computed
# field `Ex` would be required (e.g. python package nrefocus).

import numpy as np
import scattnlay
import matplotlib.pylab as plt

# weak dielectric sphere, e.g. a PMMA gel bead
n1 = 1.335
# refractive index of the surrounding medium (water)
nm = 1.333
# radius of the sphere in vacuum wavelengths
radius = 0.3
# extent of the simulation size in vacuum wavelengths
extent = 2.0
# distance where we want to have the measured field behind the sphere
# in vacuum wavelengths measured from the center of the sphere
distance = 0.5
# pixels per vacuum wavelength in the output image
resolution = 20.0 

# size parameters need to be multiplied by (2 PI nm) for the computation
twopi = 2*np.pi*nm

# There is only one sphere, no layers
x = np.ones((1, 1), dtype = np.float64)
x[0, 0] = radius*twopi

# Set the refractive index of the sphere, normalized to that of the medium
m = np.ones((1, 1), dtype = np.complex128)
m[0, 0] = n1/nm

nptsx = extent*resolution
nptsy = extent*resolution

scanx = np.linspace(-extent/2, extent/2, nptsx, endpoint=True)*twopi
scany = np.linspace(-extent/2, extent/2, nptsy, endpoint=True)*twopi

coordX, coordY = np.meshgrid(scanx, scany)
coordX.resize(nptsx*nptsy)
coordY.resize(nptsx*nptsy)
coordZ = np.ones(nptsx*nptsy, dtype=np.float64)*distance*twopi

coord = np.vstack((coordX, coordY, coordZ)).transpose()

terms, E, H = scattnlay.fieldnlay(x, m, coord)

# take the x-component of the electric field
Ex = E[:,:,0].reshape(nptsx, nptsy)

# normalize by the background field (free space propagation)
Ex /= np.exp(1j*2*np.pi*distance*nm)

# plot the phase (np.angle) of the x-component of the electric field
ax = plt.subplot(111)
mapper = plt.imshow(np.angle(Ex))
plt.colorbar(mapper, ax=ax, label="phase [rad]")
plt.title("phase retardation introduced by a dielectric sphere")
plt.show()
