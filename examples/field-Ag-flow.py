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
# E-k plane, for an spherical Ag nanoparticle.

import scattnlay
from scattnlay import fieldnlay
from scattnlay import scattnlay
from fieldplot import fieldplot

import numpy as np
import cmath
# # a)
#WL=400 #nm
#core_r = WL/20.0
#epsilon_Ag = -2.0 + 10.0j

# # b)
#WL=400 #nm
#core_r = WL/20.0
#epsilon_Ag = -2.0 + 1.0j

# c)
WL=354 #nm
core_r = WL/20.0
epsilon_Ag = -2.0 + 0.28j

# d)
#WL=367 #nm
#core_r = WL/20.0
#epsilon_Ag = -2.71 + 0.25j

index_Ag = np.sqrt(epsilon_Ag)
# n1 = 1.53413
# n2 = 0.565838 + 7.23262j
nm = 1.0

x = np.ones((2), dtype = np.float64)
x[0] = 2.0*np.pi*core_r/WL/4.0*3.0
x[1] = 2.0*np.pi*core_r/WL

m = np.ones((2), dtype = np.complex128)
m[0] = index_Ag/nm
m[1] = index_Ag/nm

print "x =", x
print "m =", m

comment='bulk-Ag-flow'
WL_units='nm'
npts = 501
factor=2.1
flow_total = 39
#flow_total = 21
#flow_total = 0
crossplane='XZ'
#crossplane='YZ'
#crossplane='XY'

# Options to plot: Eabs, Habs, Pabs, angleEx, angleHy
field_to_plot='Pabs'
#field_to_plot='angleEx'

fieldplot(x,m, WL, comment, WL_units, crossplane, field_to_plot, npts, factor, flow_total, is_flow_extend=False)

