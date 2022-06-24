#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
from matplotlib import pyplot as plt
import cmath
from scattnlay import mesomie, mie
from optical_constants import read_refractive_index_from_yaml as read_nk
import numpy as np
import scipy.io

from_disk = np.loadtxt('silver-d_perp_interpolated.txt')
omega_star_ratio = from_disk[0, :]
d_perp = from_disk[1, :] + 1j*from_disk[2, :]
from_disk = np.loadtxt('silver-d_parl_interpolated.txt')
d_parl = from_disk[1, :] + 1j*from_disk[2, :]

c = 299792458  # m/s
h_reduced = 6.5821e-16  # eV s
omega_p = 9.02  # eV
omega_p_star = 3.81  # eV
gamma = 0.022  # eV
eps_d = 1

R = 2.5
y_min = 0
y_max = 2

min_lim_omega_star_ratio = 0.87
max_lim_omega_star_ratio = 0.99

# min_lim_omega_ratio = min_lim_omega_star_ratio * omega_p_star/omega_p
# max_lim_omega_ratio = max_lim_omega_star_ratio * omega_p_star/omega_p

# 2 pi / lambda = (omega/c) /h_reduced


WL = 2*np.pi/((omega_star_ratio * omega_p_star/c)/h_reduced)*1e6  # mkm
min_WL_available = 0.1879
max_WL_available = 1.9370
WL[WL < min_WL_available] = min_WL_available
WL[WL > max_WL_available] = max_WL_available
index_Ag = read_nk('Ag-Johnson-1972.yml', WL, kind=1)
eps_Ag = index_Ag**2
print(index_Ag)


def eps_m(omega):
    return 1 - omega_p * omega_p / (omega*omega + 1j*omega*gamma)


Qext = []
Qext_mie = []
om_rat_plot = []
for i in range(len(omega_star_ratio)):
    om_star_rat = omega_star_ratio[i]
    if (om_star_rat < min_lim_omega_star_ratio
            or om_star_rat > max_lim_omega_star_ratio):
        continue
    omega = om_star_rat*omega_p_star
    WL_mkm = 2*np.pi/((omega/c)/h_reduced)*1e6
    if WL_mkm < min_WL_available or WL_mkm > max_WL_available:
        continue

    x = (omega/c) * R * 1e-9/h_reduced
    # m = cmath.sqrt(eps_m(omega))
    m = index_Ag[i, 1]
    print(x, m)
    mesomie.calc_ab(R*10,      # R in angstrem
                    x,      # xd
                    x * m,  # xm
                    1,      # eps_d
                    m * m,  # eps_m
                    d_parl[i]/2,      # d_parallel
                    d_perp[i])      # d_perp
    mesomie.calc_Q()
    mie.SetLayersSize(x)
    mie.SetLayersIndex(m)
    mie.RunMieCalculation()
    Qext.append(mesomie.GetQext())
    Qext_mie.append(mie.GetQext())
    # print(x, m, Qext[-1] - mie.GetQext())

    om_rat_plot.append(om_star_rat)

plt.plot(om_rat_plot, Qext_mie,
         label='classic', color='gray', lw=4)
plt.plot(om_rat_plot, Qext,
         label='non-classic', color='red', lw=4)
plt.legend()
# plt.yscale('log')
plt.xlim((min_lim_omega_star_ratio, max_lim_omega_star_ratio))
plt.ylim((y_min, y_max))
plt.title("R="+str(R))
plt.show()
