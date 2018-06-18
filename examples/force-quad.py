#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2018  Konstantin Ladutenko <kostyfisik@gmail.com>
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

# WIP: try to evaluate forces.

from scattnlay import fieldnlay
from scattnlay import scattnlay
import cmath
import matplotlib.pyplot as plt
import numpy as np
import scattnlay
import quadpy


index_Ag = 4.0

nm = 1.0
WL_units='nm'

x = np.ones((1), dtype = np.float64)
m = np.ones((1), dtype = np.complex128)

core_r = 120
WL = 550 

x[0] = 2.0*np.pi*core_r/WL#/4.0*3.0
m[0] = index_Ag/nm

R = x[0]*1.31

comment='bulk-NP-WL'+str(WL)+WL_units

quad_ord = 19
#quad_ord = 131

coord = quadpy.sphere.Lebedev(3).points
# coord = np.vstack((coordX, coordY, coordZ)).transpose()
def force(in_coord):
    coord = in_coord.T
    terms, Eall, Hall = fieldnlay(np.array([x]), np.array([m]), coord)
    E_all = Eall[0, :, :]
    H_all = Hall[0, :, :]
      # std::vector<double> P = (1/(2.0))
      #   *real(
      #         dot(unit,E)*vconj(E) +
      #         dot(unit,H)*vconj(H) +
      #         (-1.0/2.0)*(dot(E,vconj(E))
      #                     +dot(H,vconj(H))
      #                     )*unit
      #         );
    unit_all = coord/ R
    P = np.array([
        ( (1/(2.0))
        *np.real(
            np.dot(unit,E)*np.conj(E) +
            np.dot(unit,H)*np.conj(H) +
              (-1.0/2.0)*(np.dot(E,np.conj(E))
                          +np.dot(H,np.conj(H))
                          )*unit
              )
        )
        for unit,E,H in zip(unit_all,E_all, H_all)
        ])
    return P.T
    
# P = np.array(map(lambda n: np.linalg.norm(np.cross(Ec[n], Hc[n])).real,
#                      range(0, len(E[0]))))

val = quadpy.sphere.integrate(
    force,
    [0.0, 0.0, 0.0], R,
    quadpy.sphere.Lebedev(quad_ord)
    )
print(val)
print("Random increase of integraion sphere radius...")
val = quadpy.sphere.integrate(
    force,
    [0.0, 0.0, 0.0], R*2.718,
    quadpy.sphere.Lebedev(quad_ord)
    )
print(val)
