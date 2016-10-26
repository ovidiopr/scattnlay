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
# E-k plane, for an spherical Si-Ag-Si nanoparticle.

import scattnlay
from scattnlay import fieldnlay
from scattnlay import scattnlay
import numpy as np
import cmath
# from fieldplot import GetFlow3D
# from fieldplot import GetField
from fieldplot import fieldplot

###############################################################################
def SetXM(design):
    """ design value:
    1: AgSi - a1
    2: SiAgSi - a1, b1
    3: SiAgSi - a1, b2
    """
    epsilon_Si = 18.4631066585 + 0.6259727805j
    epsilon_Ag = -8.5014154589 + 0.7585845411j
    index_Si = np.sqrt(epsilon_Si)
    index_Ag = np.sqrt(epsilon_Ag)
    isSiAgSi=True
    isBulk = False
    if design==1:
        #36	5.62055	0	31.93	4.06	49	5.62055	500
        isSiAgSi=False
        WL=500 #nm
        core_width = 0.0 #nm Si
        inner_width = 31.93 #nm Ag
        outer_width = 4.06 #nm  Si
    elif design==2:
        #62.5	4.48866	29.44	10.33	22.73	0	4.48866	500
        WL=500 #nm
        core_width = 29.44 #nm Si
        inner_width = 10.33 #nm Ag
        outer_width = 22.73 #nm  Si
    elif design == 3:
        #81.4	3.14156	5.27	8.22	67.91	0	3.14156	500
        WL=500 #nm
        core_width = 5.27 #nm Si
        inner_width = 8.22 #nm Ag
        outer_width = 67.91 #nm  Si
    elif design==4:
        WL=800 #nm
        epsilon_Si = 13.64 + 0.047j
        epsilon_Ag = -28.05 + 1.525j
        core_width = 17.74 #nm Si
        inner_width = 23.31 #nm Ag
        outer_width = 22.95 #nm  Si
    elif design==5:
        WL=354 #nm
        core_r = WL/20.0
        epsilon_Ag = -2.0 + 0.28j   #original
        index_Ag = np.sqrt(epsilon_Ag)
        x = np.ones((1), dtype = np.float64)
        x[0] = 2.0*np.pi*core_r/WL
        m = np.ones((1), dtype = np.complex128)
        m[0] = index_Ag
        # x = np.ones((2), dtype = np.float64)
        # x[0] = 2.0*np.pi*core_r/WL/4.0*3.0
        # x[1] = 2.0*np.pi*core_r/WL
        # m = np.ones((2), dtype = np.complex128)
        # m[0] = index_Ag
        # m[1] = index_Ag
        return x, m, WL
    elif design==6:
        WL=1052 #nm
        core_r = 140.0
        #core_r = 190.0
        core_r = 204.2
        epsilon_Si = 12.7294053067+0.000835315166667j
        index_Si = np.sqrt(epsilon_Si)
        x = np.ones((1), dtype = np.float64)
        x[0] = 2.0*np.pi*core_r/WL
        m = np.ones((1), dtype = np.complex128)
        m[0] = index_Si
        # x = np.ones((2), dtype = np.float64)
        # x[0] = 2.0*np.pi*core_r/WL/4.0*3.0
        # x[1] = 2.0*np.pi*core_r/WL
        # m = np.ones((2), dtype = np.complex128)
        # m[0] = index_Ag
        # m[1] = index_Ag
        return x, m, WL


    core_r = core_width
    inner_r = core_r+inner_width
    outer_r = inner_r+outer_width

    nm = 1.0
    if isSiAgSi:
        x = np.ones((3), dtype = np.float64)
        x[0] = 2.0*np.pi*core_r/WL
        x[1] = 2.0*np.pi*inner_r/WL
        x[2] = 2.0*np.pi*outer_r/WL
        m = np.ones((3), dtype = np.complex128)
        m[0] = index_Si/nm
        m[1] = index_Ag/nm
    #    m[0, 1] = index_Si/nm
        m[2] = index_Si/nm
    else:
        # bilayer
        x = np.ones((2), dtype = np.float64)
        x[0] = 2.0*np.pi*inner_r/WL
        x[1] = 2.0*np.pi*outer_r/WL
        m = np.ones((2), dtype = np.complex128)
        m[0] = index_Ag/nm
        m[1] = index_Si/nm
    return x, m, WL


###############################################################################
#design = 1 #AgSi
#design = 2
#design = 3
# design = 4   # WL=800
# comment='SiAgSi-flow'
#design = 5   # Bulk Ag
# comment='bulk-Ag-flow'
design = 6   # WL=800
comment='Si-flow'
x, m, WL = SetXM(design)

WL_units='nm'
print "x =", x[-1]
print "m =", m
npts = 501
factor=2.1
flow_total = 39
#flow_total = 21
#flow_total = 0
#crossplane='XZ'
crossplane='XYZ'
#crossplane='YZ'
#crossplane='XY'

# Options to plot: Eabs, Habs, Pabs, angleEx, angleHy
field_to_plot='Eabs'
#field_to_plot='angleEx'
#field_to_plot='Pabs'

import matplotlib.pyplot as plt
fig, axs = plt.subplots(1,1)#, sharey=True, sharex=True)
fig.tight_layout()
fieldplot(fig, axs, x,m, WL, comment, WL_units, crossplane, field_to_plot, npts, factor, flow_total,
          subplot_label=' ',is_flow_extend=False)

fig.subplots_adjust(hspace=0.3, wspace=-0.1)

plt.savefig(comment+"-R"+str(int(round(x[-1]*WL/2.0/np.pi)))+"-"+crossplane+"-"
                    +field_to_plot+".pdf",pad_inches=0.02, bbox_inches='tight')

plt.draw()

#    plt.show()

plt.clf()
plt.close()



