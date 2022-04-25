#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2019-2021  Konstantin Ladutenko <kostyfisik@gmail.com>
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


import sys
sys.path.insert(0,'..')  # to be able to import scattnlay from the upper dir

from scattnlay import scattnlay,scattcoeffs,fieldnlay, mie

import matplotlib.pyplot as plt
import numpy as np
import cmath

from_WL = 300
to_WL = 1000
WL_points= 1000
WLs = np.linspace(from_WL, to_WL, WL_points)
index_NP = 4+0.01j

x = np.ones((1), dtype = np.float64)
m = np.ones((1), dtype = np.complex128)

core_r = 100

Qsca_vec = []
Qsca_mode_vec = []

for WL in WLs:
    x[0] = 2.0*np.pi*core_r/WL
    m[0] = index_NP
    mie.SetLayersSize(x)
    mie.SetLayersIndex(m)

    mie.SetModeNmaxAndType(-1,-1)
    mie.RunMieCalculation()
    Qsca =  mie.GetQsca()

    mie.SetModeNmaxAndType(1,0)
    mie.RunMieCalculation()
    Qsca_mode =  mie.GetQsca()

    print(np.array([Qsca]))
    Qsca_vec.append(Qsca)
    Qsca_mode_vec.append(Qsca_mode)

fig, axs2 = plt.subplots(1,1)#, sharey=True, sharex=True)
axs2.plot(WLs, Qsca_vec, color="blue")
axs2.plot(WLs, Qsca_mode_vec, color="red")
plt.savefig("spectra.pdf",pad_inches=0.02, bbox_inches='tight')
plt.show()
plt.clf()
plt.close()
