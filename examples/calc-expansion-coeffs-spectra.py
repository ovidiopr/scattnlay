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

from scattnlay import fieldnlay, scattnlay, expansioncoeffs, scattcoeffs

import numpy as np
import cmath

WL=550 #nm
core_r = 102
#core_r = 120

index_NP = 4.0
x = np.ones((1), dtype = np.float64)
m = np.ones((1), dtype = np.complex128)

import matplotlib.pyplot as plt

# dipole - 1 -- 2
# quad   - 2 -- 4
# octo   - 3 -- 8
# hex    - 4 -- 16
# 32     - 5 -- 32
npts = 151
ext = ".png"
# npts = 351
# ext = ".pdf"
x[0] = 2.0*np.pi*core_r/WL#/4.0*3.0
m[0] = index_NP
#    for mode_type, field_to_plot, WL, mode_n,  crossplane, isStream in plot_params :
comment='bulk-NP-R'+str(core_r)+'nm-WL'+str(WL)
# scattnlay function process several NP designs in one call
x = np.array([x])
m = np.array([m])
#print(x.shape)
terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(x, m)
nmax = terms
#print(nmax)
terms, an, bn = scattcoeffs(x, m, nmax)
#print(terms, an)
terms, aln, bln, cln, dln = expansioncoeffs(x, m, nmax)
print(terms, dln)
# fig, axs = plt.subplots(1,1)#, sharey=True, sharex=True)
# fig.tight_layout()
# fig.subplots_adjust(hspace=0.3, wspace=-0.1)

# # plt.savefig("Egor3/"+"%02d"%(i)+
# #             comment+"-R"+str(int(round(x[-1]*WL/2.0/np.pi)))+"-"+crossplane+"-"
# #                 +field_to_plot+"-mode"+mode+mt+st+ext,pad_inches=0.02, bbox_inches='tight')
# plt.clf()
# plt.close()

#print("end")
