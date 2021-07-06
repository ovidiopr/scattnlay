#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2021  Konstantin Ladutenko <kostyfisik@gmail.com>
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

from scattnlay import mie
import matplotlib.pyplot as plt
import numpy as np
from optical_constants import read_refractive_index_from_yaml as get_index

def gauss(x, mu, sigma):
    return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2))

from_WL = 300
to_WL = 1100
WL_points= 100
WLs = np.linspace(from_WL, to_WL, WL_points)

from_r = 40/2.
to_r =80/2.
r_points = 20
all_r = np.linspace(from_r, to_r, r_points)
r_mean = 58.3/2.
# r_mean = 50/2.
r_std = 6.3/2.
r_weights = gauss(all_r, r_mean,r_std)/len(all_r)

plt.plot(all_r, r_weights )
plt.xlabel("R, nm")
plt.ylabel("amount")

index_SiO2 = get_index("refractiveindex_info/SiO2-Gao.yml", WLs, units='nm')
# index_Au = get_index("refractiveindex_info/Au-McPeak.yml", WLs, units='nm')
index_Au = get_index("refractiveindex_info/Au-Johnson.yml", WLs, units='nm')
index_TiO2 = get_index("r"
                       "efractiveindex_info/TiO2-Sarkar.yml", WLs, units='nm')

index_SiO2 *= 0; index_SiO2 += 1.45
index_TiO2[:,1] += 0.0j
# index_Au[:,1] += 1.5j

x = np.ones((3), dtype = np.float64)
m = np.ones((3), dtype = np.complex128)

core_r = 5
inner_shell_h = 10+20
outer_shell_h = 10
host_media = 1.33

Qext_core_shell = np.zeros(len(WLs))
Qext_3l = np.zeros(len(WLs))

for i in range(len(WLs)):
    WL = WLs[i]
    for j in range(len(all_r)):
        # core_r = all_r[j]
        # weight = r_weights[j]
        weight = 1/len(r_weights)
        # print(core_r)
        x = host_media*2.0*np.pi/WL*np.array([core_r,
                                   core_r+inner_shell_h,
                                   core_r+inner_shell_h+outer_shell_h])
        m = np.array([index_SiO2[i][1], index_Au[i][1],
                      index_TiO2[i][1]]
                     )/host_media
        # print(x, m)
        mie.SetLayersSize(x)
        mie.SetLayersIndex(m)
        mie.RunMieCalculation()
        Qext_3l[i] += mie.GetQext()*weight

        x = host_media*2.0*np.pi/WL*np.array([core_r,
                          core_r+inner_shell_h])
        m = np.array([index_SiO2[i][1], index_Au[i][1]])/host_media
        mie.SetLayersSize(x)
        mie.SetLayersIndex(m)
        mie.RunMieCalculation()
        Qext_core_shell[i] += mie.GetQext()*weight

comsol_spectra = np.array([[0.420000000000000,2.35836000000000e-15],
                  [0.440000000000000,2.27000000000000e-15],
                  [0.460000000000000,2.21146900000000e-15],
                  [0.480000000000000,2.21744500000000e-15],
                  [0.500000000000000,2.49989500000000e-15],
                  [0.520000000000000,3.36257000000000e-15],
                  [0.540000000000000,3.88983000000000e-15],
                  [0.560000000000000,4.03982000000000e-15],
                  [0.580000000000000,3.23889000000000e-15],
                  [0.600000000000000,3.01499000000000e-15],
                  [0.620000000000000,2.13147000000000e-15],
                  [0.640000000000000,9.02930000000000e-16],
                  [0.660000000000000,4.49688000000000e-16],
                  [0.680000000000000,2.93514000000000e-16],
                  [0.700000000000000,2.19381000000000e-16],
                  [0.720000000000000,1.85272000000000e-16],
                  [0.740000000000000,1.74517000000000e-16],
                  [0.760000000000000,1.54702000000000e-16],
                  [0.780000000000000,1.51191000000000e-16],
                  [0.800000000000000,1.58785200000000e-16],
                  [0.820000000000000,1.74967600000000e-16 ]  ])

fig, axs2 = plt.subplots(1,1)#, sharey=True, sharex=True)
# axs2.plot(WLs, Qext_3l, color="purple")
axs2.plot(WLs, Qext_core_shell, color="lime", label="Mie, layered")
axs2.plot(comsol_spectra[:,0]*1000, comsol_spectra[:,1]/np.pi/(35e-9**2)*3, color="black", label="Comsol, \nSiO2 with Au NP coating")
axs2.legend()
axs2.set_xlabel("WL, nm")
axs2.set_ylabel("Extinction, a.u.")
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
