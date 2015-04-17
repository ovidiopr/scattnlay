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

# This test case calculates the electric field along three 
# points, for an spherical silver nanoparticle embedded in glass.

# Refractive index values correspond to a wavelength of
# 400 nm. Maximum of the surface plasmon resonance (and,
# hence, of electric field) is expected under those
# conditions.

from scattnlay import fieldnlay
import numpy as np

n1 = 1.53413
n2 = 0.565838 + 7.23262j
nm = 1.3205

x = np.ones((1, 3), dtype = np.float64)
x[0, 0] = 2.0*np.pi*nm*0.05/1.064
x[0, 1] = 2.0*np.pi*nm*0.06/1.064
x[0, 2] = 2.0*np.pi*nm*0.07/1.064

m = np.ones((1, 3), dtype = np.complex128)
m[0, 0] = n1/nm
m[0, 1] = n2/nm
m[0, 2] = 1.0

coord = np.zeros((3, 3), dtype = np.float64)
coord[0, 0] = x[0, 0]/2.0
coord[1, 0] = (x[0, 0] + x[0, 1])/2.0
coord[2, 0] = 1.5*x[0, 1]

terms, E, H = fieldnlay(x, m, coord)

Er = np.absolute(E)

# |E|/|Eo|
Eh = np.sqrt(Er[0, :, 0]**2 + Er[0, :, 1]**2 + Er[0, :, 2]**2)

print "x =", x
print "m =", m
print np.vstack((coord[:, 0], Eh)).transpose()

