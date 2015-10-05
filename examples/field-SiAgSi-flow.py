#!/usr/bin/env python
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
#    using it, cite the following reference:
#    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by
#        a multilayered sphere," Computer Physics Communications,
#        vol. 180, Nov. 2009, pp. 2348-2354.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# This test case calculates the electric field in the 
# E-k plane, for an spherical Si-Ag-Si nanoparticle. Core radius is 17.74 nm,
# inner layer 23.31nm, outer layer 22.95nm. Working wavelength is 800nm, we use
# silicon epsilon=13.64+i0.047, silver epsilon= -28.05+i1.525

import scattnlay
from scattnlay import fieldnlay
from scattnlay import scattnlay
from fieldplot import fieldplot
import numpy as np
import cmath


epsilon_Si = 13.64 + 0.047j
epsilon_Ag = -28.05 + 1.525j

# epsilon_Si = 2.0 + 0.047j
# epsilon_Ag = -2.0 + 1.525j

# air = 1
# epsilon_Si = air*2
# epsilon_Ag = air*2


index_Si = np.sqrt(epsilon_Si)
index_Ag = np.sqrt(epsilon_Ag)

print(index_Si)
print(index_Ag)
# # Values for 800 nm, taken from http://refractiveindex.info/
# index_Si = 3.69410 + 0.0065435j
# index_Ag = 0.18599 + 4.9886j

WL=800 #nm
core_width = 17.74 #nm Si
inner_width = 23.31 #nm Ag
outer_width = 22.95 #nm  Si

core_r = core_width
inner_r = core_r+inner_width
outer_r = inner_r+outer_width

# n1 = 1.53413
# n2 = 0.565838 + 7.23262j
nm = 1.0

x = np.ones((3), dtype = np.float64)
x[0] = 2.0*np.pi*core_r/WL
x[1] = 2.0*np.pi*inner_r/WL
x[2] = 2.0*np.pi*outer_r/WL

m = np.ones((3), dtype = np.complex128)
m[0] = index_Si/nm
m[1] = index_Ag/nm
m[2] = index_Si/nm

print "x =", x
print "m =", m

npts = 501
factor=2.2
flow_total = 21


crossplane='XZ'
#crossplane='YZ'
#crossplane='XY'

# Options to plot: Eabs, Habs, Pabs, angleEx, angleHy
field_to_plot='Eabs'
#field_to_plot='angleEx'
comment='SiAgSi-absorber-flow'
WL_units='nm'

import matplotlib.pyplot as plt
fig, axs = plt.subplots(1,1)#, sharey=True, sharex=True)
fig.tight_layout()
fieldplot(fig, axs, x,m, WL, comment, WL_units, crossplane, field_to_plot, npts, factor, flow_total,
          subplot_label=' ',is_flow_extend=False, outline_width=1.5)

#fieldplot(x,m, WL, comment, WL_units, crossplane, field_to_plot, npts, factor, flow_total, is_flow_extend=False)

# for ax in axs:
#     ax.locator_params(axis='x',nbins=5)
#     ax.locator_params(axis='y',nbins=5)

fig.subplots_adjust(hspace=0.3, wspace=-0.1)

plt.savefig(comment+"-R"+str(int(round(x[-1]*WL/2.0/np.pi)))+"-"+crossplane+"-"
                    +field_to_plot+".pdf",pad_inches=0.02, bbox_inches='tight')

plt.draw()

#    plt.show()

plt.clf()
plt.close()


