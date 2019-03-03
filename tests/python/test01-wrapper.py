#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2017 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#
#    This file is part of scattnlay
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
#    using it, cite at least one of the following references:
#    [1] O. Peña and U. Pal, "Scattering of electromagnetic radiation by
#        a multilayered sphere," Computer Physics Communications,
#        vol. 180, Nov. 2009, pp. 2348-2354.
#    [2] K. Ladutenko, U. Pal, A. Rivera, and O. Peña-Rodríguez, "Mie
#        calculation of electromagnetic near-field for a multilayered
#        sphere," Computer Physics Communications, vol. 214, May 2017,
#        pp. 225-230.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# This is a test against the program n-mie (version 3a) for the test case
# distributed by them (extended for x up to 100)
# n-mie is based in the algorithm described in:
# Wu Z.P., Wang Y.P.
# Electromagnetic scattering for multilayered spheres:
# recursive algorithms
# Radio Science 1991. V. 26. P. 1393-1401.
# Voshchinnikov N.V., Mathis J.S.
# Calculating Cross Sections of Composite Interstellar Grains
# Astrophys. J. 1999. V. 526. #1. 

# The test consist in 5 layers with the following parameters
# m1=1.8 i1.7
# m2=0.8 i0.7
# m3=1.2 i0.09
# m4=2.8 i0.2
# m5=1.5 i0.4

# v1/Vt=0.1
# v2/Vt=0.26
# v3/Vt=0.044
# v4/Vt=0.3666
import scattnlay

# import os
# path = os.path.dirname(scattnlay.__file__)
print(scattnlay.__file__)

#from scattnlay import scattnlay
#import scattnlay
import numpy as np
# import os
# import inspect
# inspect.getfile(scattnlay)

x = np.ones((400, 5), dtype = np.float64)
x[:, 4] = np.arange(0.25, 100.25, 0.25)
x[:, 0] = 0.1**(1.0/3.0)*x[:, 4]
x[:, 1] = 0.36**(1.0/3.0)*x[:, 4]
x[:, 2] = 0.404**(1.0/3.0)*x[:, 4]
x[:, 3] = 0.7706**(1.0/3.0)*x[:, 4]

m = np.ones((400, 5), dtype = np.complex128)
m[:, 0] *= 1.8 + 1.7j
m[:, 1] *= 0.8 + 0.7j
m[:, 2] *= 1.2 + 0.09j
m[:, 3] *= 2.8 + 0.2j
m[:, 4] *= 1.5 + 0.4j

terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay.scattnlay_wrapper(x, m)

result = np.vstack((x[:, 4], Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo)).transpose()

try:
    import matplotlib.pyplot as plt

    plt.figure(1)
    plt.subplot(311)
    plt.plot(x[:, 4], Qext, 'k')
    plt.ylabel('Qext')

    plt.subplot(312)
    plt.plot(x[:, 4], Qsca, 'r')
    plt.ylabel('Qsca')

    plt.subplot(313)
    plt.plot(x[:, 4], Albedo, 'g')
    plt.ylabel('Albedo')

    plt.xlabel('X')
    
    #plt.show()
finally:
    np.savetxt("test01.txt", result, fmt = "%.5f")
    print result

