#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2015 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#    Copyright (C) 2013-2015  Konstantin Ladutenko <kostyfisik@gmail.com>
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
#    using it, cite at least one of the following references:
#    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by
#        a multilayered sphere," Computer Physics Communications,
#        vol. 180, Nov. 2009, pp. 2348-2354.
#    [2] K. Ladutenko, U. Pal, A. Rivera, and O. Pena-Rodriguez, "Mie
#        calculation of electromagnetic near-field for a multilayered
#        sphere," Computer Physics Communications, vol. 214, May 2017,
#        pp. 225-230.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http:#www.gnu.org/licenses/>.

import scattnlay
from scattnlay import fieldnlay
from scattnlay import scattnlay
import numpy as np
from matplotlib import pyplot as plt

npts = 151
factor = 3. # plot extent compared to sphere radius
index_H2O = 1.33+0.j

WL = 0.532 #mkm
total_r = 1 #mkm
isMP = False
# isMP = True
nm = 1.0 # pihost medium
x = 2.0 * np.pi * np.array([total_r / 4.0 * 3.0, total_r], dtype=np.float64) / WL
m = np.array((index_H2O, index_H2O), dtype=np.complex128) / nm

print("x =", x)
print("m =", m)
terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(
    np.array([x]), np.array([m]))
print("Qsca = " + str(Qsca))
terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(
    np.array([x]), np.array([m]), mp=True)
print("mp Qsca = " + str(Qsca))

scan = np.linspace(-factor*x[-1], factor*x[-1], npts)
zero = np.zeros(npts*npts, dtype = np.float64)

coordX, coordZ = np.meshgrid(scan, scan)
coordX.resize(npts * npts)
coordZ.resize(npts * npts)
coordY = zero

terms, E, H = fieldnlay(
    np.array([x]), np.array([m]),
    coordX, coordY, coordZ,
    mp=isMP
)
Ec = E[0, :, :]
Er = np.absolute(Ec)
Eabs2 = (Er[:, 0]**2 + Er[:, 1]**2 + Er[:, 2]**2)
Eabs_data = np.resize(Eabs2, (npts, npts))
label = r'$|E|^2$'
plt.imshow(Eabs_data,
           cmap='jet',
           vmin=0., vmax=14

           )
print(np.min(Eabs_data), np.max(Eabs_data))
mp = ''
if isMP: mp = '_mp'
plt.savefig("R"+str(total_r)+"mkm"+mp+".jpg")
# plt.show()
