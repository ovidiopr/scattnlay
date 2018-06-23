#!/usr/bin/env python
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

R_st = x[0]*2
#R_st = 0.31

dx = R_st*4.0
charge = 1.0

comment='bulk-NP-WL'+str(WL)+WL_units

#quad_ord = 3
# quad_ord = 19
quad_ord = 131

#coord = quadpy.sphere.Lebedev(3).points
# coord = np.vstack((coordX, coordY, coordZ)).transpose()
def field(coord):
    E = []
    for rr in coord:
        shift = np.array([dx, 0.0, 0.0])
        unit = (rr-shift)/np.linalg.norm(rr-shift)
        norm = np.linalg.norm(rr-shift)
        amp = charge/(4*np.pi*(norm**2))
        Eloc = amp*unit
        
        shift = np.array([0.0, 0.0, 0.0])
        unit = (rr-shift)/np.linalg.norm(rr-shift)
        norm = np.linalg.norm(rr-shift)
        amp = charge/(4*np.pi*(norm**2))
        Eloc += amp*unit
        
        E.append(Eloc)
    E = np.array(E)
    return E

def gauss(in_coord):
    coord = in_coord.T
    E_all = field(coord)
    unit_all = coord/ R
    P = np.array([
        np.dot(E,unit)
        for unit,E in zip(unit_all,E_all)
        ])
    return P.T

def dipole(coord):
    H = np.array([0.0,0.0,0.0]*len(coord))
    E = field(coord)
    return E,H

def force(in_coord):
    coord = in_coord.T
    terms, Eall, Hall = fieldnlay(np.array([x]), np.array([m]), coord)#, mode_n=-1, mode_type=0)
    E_all = Eall[0, :, :]
    H_all = Hall[0, :, :]

    E_all, H_all = dipole(coord)
    # print(coord)
    # print(E_all)
    # print(H_all)

    unit_all = coord/ R
    P = np.array([
        ( (1.0/(2.0)*(4.0*np.pi))
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

def poynting(in_coord):
    coord = in_coord.T
    terms, Eall, Hall = fieldnlay(np.array([x]), np.array([m]), coord)
    E_all = Eall[0, :, :]
    H_all = Hall[0, :, :]
    unit_all = coord/ R
    P = np.array([
        ( ( 1.0/(2.0) )
            *np.real(
                np.cross(E,np.conj(H))
            )
        )
        for unit,E,H in zip(unit_all,E_all, H_all)
        ])
    return P.T


# P = np.array(map(lambda n: np.linalg.norm(np.cross(Ec[n], Hc[n])).real,
#                      range(0, len(E[0]))))
R=R_st
val = quadpy.sphere.integrate(
    force
#    poynting
    ,
    [0.0, 0.0, 0.0], R,
    quadpy.sphere.Lebedev(quad_ord)
    )
print(val)
print("Random increase of integraion sphere radius...")
R=R_st*3.0
val = quadpy.sphere.integrate(
    force
   # poynting
    ,
    [0.0, 0.0, 0.0], R,
    quadpy.sphere.Lebedev(quad_ord)
    )
print(val)

R=R_st
print("\n\nCheck Gauss law")
val = quadpy.sphere.integrate(gauss,
    [0.0, 0.0, 0.0], R,
    quadpy.sphere.Lebedev(quad_ord)
    )
print("Charge: ",val)
print("Random increase of integraion sphere radius...")
R=R_st*3
val = quadpy.sphere.integrate(gauss,
        [0.0, 0.0, 0.0], R,
    quadpy.sphere.Lebedev(quad_ord)
    )
print("Charge: ",val)
