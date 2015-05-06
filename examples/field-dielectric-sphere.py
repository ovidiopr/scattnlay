#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2015 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
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
# small dielectric sphere.
import scattnlay
import os
from scattnlay import fieldnlay
import numpy as np

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def is_test_coord_passed(x,m,coord):
    terms, E, H = fieldnlay(x, m, coord)
    Er = np.absolute(E)
    Eabs = np.sqrt(Er[0, :, 0]**2 + Er[0, :, 1]**2 + Er[0, :, 2]**2)
    analytic_E = (3/(m[0,0]**2+2)).real
    for value in Eabs:
        #print(value, value-analytic_E)
        if ( value-analytic_E > 10e-7 ):
            print(bcolors.FAIL+"Test failed: value="+str(value)+" for m="+str(m[0,0])
                   +" instead of analytic Eabs="+str(analytic_E))
            print("Coords",coord)
            return False
    return True

def is_test_all_coord_passed(x,m):
    npts = 5
    scan = np.linspace(0.999*x[0, 0], -0.999*x[0, 0], npts)
    zero = np.zeros(npts, dtype = np.float64)
    coordZ = np.vstack((zero, zero, scan)).transpose()
    coordY = np.vstack((zero, scan, zero)).transpose()
    coordX = np.vstack((scan, zero, zero)).transpose()
    if (is_test_coord_passed(x,m,coordX)
        and is_test_coord_passed(x,m,coordY)
        and is_test_coord_passed(x,m,coordY)):
        return True
    return False

def test_sphere():
    path = os.path.dirname(scattnlay.__file__)
    print(bcolors.HEADER+"===== Small dielectric sphere test ====="+bcolors.ENDC)
    print("Test for python module of Scattnlay: "+scattnlay.__file__)
    #Set the sphere
    x = np.ones((1, 1), dtype = np.float64)
    x[0, 0] = 0.0001
    m = np.ones((1, 1), dtype = np.complex128)
    m[0, 0] = 1.0
    delta_m = 0.712
    for n in xrange(0,15):
        m[0,0] = 1.0 + delta_m*n
        print("Testing index m="+str(m[0,0]))
        if not is_test_all_coord_passed(x,m):
            print(bcolors.FAIL+"Test for dielectric sphere failed!"+bcolors.ENDC)
            return False
    print(bcolors.OKGREEN+"All tests for dielectric sphere passed!"+bcolors.ENDC)
    return True

test_sphere()



