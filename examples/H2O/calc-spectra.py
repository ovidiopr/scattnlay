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
sys.path.insert(0,'../..')  # to be able to import scattnlay from the upper dir

from scattnlay import scattnlay,scattcoeffs,fieldnlay

import matplotlib.pyplot as plt
import numpy as np
import cmath

from optical_constants import read_refractive_index_from_yaml as get_index

from_rWL = 0.01
to_rWL = 5  # limit from H2O-Hale.yml data
step_rWL = 0.01
rWLs = np.arange(from_rWL, to_rWL+step_rWL/2., step_rWL);
WLs = 1/rWLs #mkm

index_H2O = get_index('H2O-Hale.yml', WLs, "mkm")

print(index_H2O)

x = np.ones((1), dtype = np.float64)
m = np.ones((1), dtype = np.complex128)


core_r = 1 #mkm


Qext_vec = []

for i in range(len(WLs)):
    WL = WLs[i]
    x[0] = 2.0*np.pi*core_r/WL#/4.0*3.0
    m[0] = index_H2O[:,1][i]
    
    terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(
        np.array(x), np.array(m),
        mp=True
        # mp=False
    )
    print(np.array([Qext]))
    Qext_vec.append(Qext)

fig, axs = plt.subplots(1,1)#, sharey=True, sharex=True)
axs.plot(rWLs, Qext_vec, color="black")
plt.ylim(0, 4.3) 
axs.set_xlabel("$1/\lambda, \mu m^{-1}$")
axs.set_ylabel("$Q_{ext}$")
plt.title("Scattnlay")
plt.savefig("spectra.pdf",pad_inches=0.02, bbox_inches='tight')
plt.show()
plt.clf()
plt.close()
