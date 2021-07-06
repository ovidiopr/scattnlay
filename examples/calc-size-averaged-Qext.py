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

fill_factor = 0.8

from_WL = 300
to_WL = 1100
WL_points= 100
WLs = np.linspace(from_WL, to_WL, WL_points)

r_points = 200

from_r = 40/2.
to_r =80/2.
all_r = np.linspace(from_r, to_r, r_points)
r_mean = 58.3/2.
r_std = 6.3/2.


# from_r = 1.
# to_r =15.
# all_r = np.linspace(from_r, to_r, r_points)
# r_mean = 6
# r_std = 2

core_r = 30
inner_shell_h = 8
outer_shell_h = 10
host_media = 1.33

r_weights = gauss(all_r, r_mean,r_std)/len(all_r)

plt.plot(all_r, r_weights )
plt.xlabel("R, nm")
plt.ylabel("amount")

index_SiO2 = get_index("refractiveindex_info/SiO2-Gao.yml", WLs, units='nm')
# index_Au = get_index("refractiveindex_info/Au-McPeak.yml", WLs, units='nm')
index_Au = get_index("refractiveindex_info/Au-Johnson.yml", WLs, units='nm')
index_TiO2 = get_index("refractiveindex_info/TiO2-Sarkar.yml", WLs, units='nm')

# index_TiO2[:,1] += 0.5j
# index_Au[:,1] = index_Au[:,1]* fill_factor + (1-fill_factor)

print("Au index before correction, max = ", np.max(np.imag(index_Au[:,1])))
# Taking into account gold free electrons damping
# contributed by surface scattering and bulk dumping
# See eq 1 in [1] -> doi: https://doi.org/10.1186/s11671-018-2670-7
eps_exp = index_Au[:,1]**2
c=299792458  # speed of light
h= 4.135667516e-15  # eV*s, Planck constant
w = h*c/(WLs*1e-9)  # eV, frequency
w_p =8.55  # eV, gold plasmon frequency
g_b = 18.4e-3  # eV, bulk dumping
v_f = 1.4e6  # m/s, Fermi velocity of electrons in gold
A = 1.33  # fit parameter for 16 nm gold shell, see Table 2 in [1]
L_b = (4.*((core_r+inner_shell_h)**3 - core_r**3)/
       (3.*((core_r+inner_shell_h)**2 + core_r**2)))
# g_s = v_f/L_b  # eq 2 in [1]
g_s = h*A*v_f/(inner_shell_h*1e-9)  # eq 4 in [1]
# g_b *= 0.2
# g_s *= 0.5
eps_Au = (eps_exp
          +
          w_p**2 / ( w * (w + 1j*g_b) )
          -
          w_p**2 / ( w * (w + 1j*(g_b+g_s)) )
          )
# index_Au[:,1] = np.sqrt(eps_Au)
index_Au[:,1] += 1.6j
print(f"L_b={L_b}")
# print(w)
print(f"g_s={g_s}, g_b={g_b} ")
print("Au index after, max = ", np.max(np.imag(index_Au[:,1])))
x = np.ones((3), dtype = np.float64)
m = np.ones((3), dtype = np.complex128)


Qext_core_shell = np.zeros(len(WLs))
Qext_3l = np.zeros(len(WLs))

for i in range(len(WLs)):
    WL = WLs[i]
    for j in range(len(all_r)):
        core_r = all_r[j]
        # inner_shell_h = all_r[j]
        weight = r_weights[j]
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


fig, axs2 = plt.subplots(1,1)#, sharey=True, sharex=True)
axs2.plot(WLs, Qext_3l, color="purple")
axs2.plot(WLs, Qext_core_shell, color="lime")
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
