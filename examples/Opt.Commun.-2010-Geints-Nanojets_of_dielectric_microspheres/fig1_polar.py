#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2021 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#    Copyright (C) 2013-2021  Konstantin Ladutenko <kostyfisik@gmail.com>
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
#    using it, cite at least one of the following references:
#    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by
#        a multilayered sphere," Computer Physics Communications,
#        vol. 180, Nov. 2009, pp. 2348-2354.
#    [2] K. Ladutenko, U. Pal, A. Rivera, and O. Pena-Rodriguez, "Mie
#        calculation of electromagnetic near-field for a multilayered
#        sphere," Computer Physics Communications, vol. 214, May 2017,
#        pp. 225-230.
#
import numpy as np
import matplotlib.pyplot as plt
from scattnlay import mie, mie_mp

# npts = 151/2
npts = 51/2
factor = 2  # plot extent compared to sphere radius
total_r = 5  # mkm
isMP = False
# isMP = True

terms_in = 210
# terms = -1


from_theta = 0
to_theta = np.pi*2
outer_arc_points = int(abs(to_theta-from_theta)*npts)
# outer_arc_points = 600


index_H2O = 1.33+(1e-6)*1j
# index_H2O = 1.001

WL = 0.532 #mkm

mp = ''
if isMP: mp = '_mp'


nm = 1.0  # host medium
x = 2.0 * np.pi * np.array([total_r/2, total_r], dtype=np.float64) / WL
m = np.array((index_H2O, index_H2O), dtype=np.complex128) / nm

# x = 2.0 * np.pi * np.array([total_r], dtype=np.float64) / WL
# m = np.array((index_H2O), dtype=np.complex128) / nm

from_r = 0.01*x[-1]
to_r = x[-1]*factor
r_points = int(outer_arc_points/abs(to_theta-from_theta))

# nmax = int(np.max(x*factor + 11 * (x*factor)**(1.0 / 3.0) + 1))
nmax = -1

print("x =", x)
print("m =", m)
mie.SetLayersSize(x)
mie.SetLayersIndex(m)
mie.RunMieCalculation()
Qsca =  mie.GetQsca()
terms = mie.GetMaxTerms()
print("   Qsca = " + str(Qsca)+" terms = "+str(terms))
mie_mp.SetLayersSize(x)
mie_mp.SetLayersIndex(m)
mie_mp.RunMieCalculation()
Qsca =  mie_mp.GetQsca()
terms = mie_mp.GetMaxTerms()

print("mp Qsca = " + str(Qsca)+" terms = "+str(terms))
# exit(1)

# mie.SetFieldCoords(xp, yp, zp)
# mie.RunFieldCalculation()
# terms = mie.GetMaxTerms()
# E = mie.GetFieldE()
# H = mie.GetFieldH()

theta_all = np.linspace(from_theta, to_theta, outer_arc_points);
r_all = np.linspace(from_r, to_r, r_points)

theta = []
r = []
for i in range(len(r_all)):
    for j in range(len(theta_all)):
        theta.append(theta_all[j])
        r.append(r_all[i])
if isMP:
    mie_mp.RunFieldCalculationPolar(outer_arc_points, r_points, from_r, to_r, from_theta, to_theta, 0, 0, True, terms_in)
    Eabs = (mie_mp.GetFieldEabs())**2
    terms = mie_mp.GetMaxTerms()
else:
    mie.RunFieldCalculationPolar(outer_arc_points, r_points, from_r, to_r, from_theta, to_theta, 0, 0, True, terms_in)
    Eabs = (mie.GetFieldEabs())**2
    terms = mie.GetMaxTerms()
print(mp, "min(Eabs)=", np.min(Eabs)," max(Eabs)=", np.max(Eabs)," terms = "+str(terms), ' size=', Eabs.size)


fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))

vmax = 14
Eabs[Eabs>vmax] = vmax
pos = ax.tricontourf(theta, r, Eabs,
                     levels=1000,
               cmap='gnuplot',
               )
plt.colorbar(pos)
ax.yaxis.grid(False)
ax.xaxis.grid(False)

# ax.plot(theta,r, 'ro', ms=0.1)
plt.savefig("R"+str(total_r)+"mkm"+mp+"_polar.jpg",
            dpi=600
            )
# plt.show()
