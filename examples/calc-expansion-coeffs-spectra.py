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

# This script reproduces Fig.2a from "Nonradiating anapole modes in
# dielectric nanoparticles" by Miroshnichenko, Andrey E. et al in  
# Nature Communications DOI:10.1038/ncomms9069


from scattnlay import fieldnlay, scattnlay, expansioncoeffs, scattcoeffs
import cmath
import matplotlib.pyplot as plt
import numpy as np

x = np.ones((1), dtype = np.float64)
m = np.ones((1), dtype = np.complex128)

WL=550 #nm
core_r = 180
index_NP = 4.0

from_R = 120/2.0
to_R = 240/2.0

npts = 151
ext = ".png"
npts = 351
# ext = ".pdf"
comment='bulk-NP-WL'+str(WL)

val_all = []

all_R = np.linspace(from_R, to_R, npts)
for core_r in all_R:
    x[0] = 2.0*np.pi*core_r/WL
    m[0] = index_NP
    terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(np.array([x]), np.array([m]))
    terms, an, bn = scattcoeffs(np.array([x]), np.array([m]),terms)

    terms, aln, bln, cln, dln = expansioncoeffs(np.array([x]), np.array([m]), terms)
    order = 0
    val_all.append([#Qsca,
                    np.abs(aln[0][-1][order]),np.abs(dln[0][0][order])
        #,np.abs(an[0][order])
    ])
#print( )
val_all = np.array(val_all)
#print()
#print(terms, np.abs(dln))
fig, axs = plt.subplots(1,1)#, sharey=True, sharex=True)
fig.tight_layout()
fig.subplots_adjust(hspace=0.3, wspace=-0.1)
all_R = all_R*2
plt.plot(all_R, val_all[:,0],label="electric dipole, scatt.")
plt.plot(all_R, val_all[:,1],label="electric dipole, internal")
plt.legend()
#plt.plot(all_R, val_all[:,2])
plt.ylim(0,4)
plt.xlabel("D, nm")
plt.ylabel("Mie coefficients")
plt.savefig(comment+"-R"+str(int(round(x[-1]*WL/2.0/np.pi)))+ext,
            pad_inches=0.02, bbox_inches='tight')
plt.show()
plt.clf()
plt.close()

#print("end")
