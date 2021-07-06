#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2019  Konstantin Ladutenko <kostyfisik@gmail.com>
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

from scattnlay import scattnlay,scattcoeffs,fieldnlay

import matplotlib.pyplot as plt
import numpy as np
import cmath


from_WL = 400
to_WL = 800
WL_points= 100
WLs = np.linspace(from_WL, to_WL, WL_points)
index_NP = 1.5+0.3j

x = np.ones((1), dtype = np.float64)
m = np.ones((1), dtype = np.complex128)


core_r = 100


Qsca_vec = []
core_r_vec = []
an_vec = []
bn_vec = []

for WL in WLs:
    x[0] = 2.0*np.pi*core_r/WL#/4.0*3.0
    m[0] = index_NP
    terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(
        np.array(x), np.array(m),
        mp=True
    )
    print(np.array([Qsca]))
    terms, an, bn = scattcoeffs(x, m,24)
    # Qsca_vec.append(Qsca*np.pi*core_r**2*1e-5)
    Qsca_vec.append(Qsca)#*np.pi*core_r**2*1e-5)
    core_r_vec.append(core_r)
    an_vec.append(np.abs(an)[0])
    bn_vec.append(np.abs(bn)[0])

an_vec = np.array(an_vec)
bn_vec = np.array(bn_vec)
# print(an_vec)
fig, axs2 = plt.subplots(1,1)#, sharey=True, sharex=True)
axs2.plot(WLs, Qsca_vec, color="black")
# axs.set_xlabel("D, nm")
# axs.set_ylabel("$Q_{sca}$")
# axs2 = axs.twinx()
# axs2.plot(np.array(core_r_vec)*2,an_vec[:,0],"b.",lw=0.8, markersize=1.9,label="$a_0$")
# axs2.plot(np.array(core_r_vec)*2,bn_vec[:,0],"b-", markersize=1.9,label="$b_0$")
# axs2.plot(np.array(core_r_vec)*2,an_vec[:,1],"g.",lw=0.8, markersize=1.9,label="$a_1$")
# axs2.plot(np.array(core_r_vec)*2,bn_vec[:,1],"g-", markersize=1.9,label="$b_1$")
# axs2.legend(loc="upper right")
# axs2.tick_params('y', colors='black')
# axs2.set_ylim(0,1)
# axs2.set_ylabel("Mie",color="black")
plt.savefig("spectra.pdf",pad_inches=0.02, bbox_inches='tight')
plt.show()
plt.clf()
plt.close()
