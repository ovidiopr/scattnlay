#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import cmath

import numpy as np
import scipy.io
from matplotlib import pyplot as plt
from optical_constants import read_refractive_index_from_yaml as read_nk

from scattnlay import mesomie, mie
# shell 0.4-0.6 nm
# omega_p goes down
shell_h = 0.4

from_disk = np.loadtxt('silver-d_perp_interpolated.txt')
omega_star_ratio = from_disk[0, :]
d_perp = from_disk[1, :] + 1j*from_disk[2, :]
from_disk = np.loadtxt('silver-d_parl_interpolated.txt')
d_parl = from_disk[1, :] + 1j*from_disk[2, :]

c = 299792458  # m/s
h_reduced = 6.5821e-16  # eV s
omega_p = 9.02  # eV
omega_p_star = 3.81  # eV
gamma = 0.15  # eV
eps_inf_drud = 4.65
# eps_inf = 4.65

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
# print(index_Ag)

factor = 1


def eps_m(omega, eps_inf, omega_p_local):
    return eps_inf - omega_p_local * omega_p_local / (omega*omega + 1j*omega*gamma*factor)


def eps_inf(omega, eps_exp):
    return eps_exp + omega_p * omega_p / (omega*omega + 1j*omega*gamma*factor)


Qext = []
Qext_mie = []
Qext_drude_nc = []
om_rat_plot = []
eps_inf_drude = []
eps_m_drude = []
for i in range(len(omega_star_ratio)):
    om_star_rat = omega_star_ratio[i]
    if (om_star_rat < min_lim_omega_star_ratio
            or om_star_rat > max_lim_omega_star_ratio):
        continue
    omega = om_star_rat*omega_p_star
    WL_mkm = 2*np.pi/((omega/c)/h_reduced)*1e6
    if WL_mkm < min_WL_available or WL_mkm > max_WL_available:
        continue

    x_const = (omega/c) * 1e-9/h_reduced
    x = R * x_const
    m = index_Ag[i, 1]
    eps_m_drude.append(m**2)
    eps_inf_drude.append(eps_inf(omega, m**2))
    m_drude = cmath.sqrt(eps_m(omega, eps_inf(omega, m**2), omega_p))
    # m_drude = cmath.sqrt(eps_m(omega, eps_inf_drud, omega_p))
    # m_drude = cmath.sqrt(eps_inf(omega, m**2))
    print(x, m)
    # m = m_drude
    mesomie.calc_ab(R*10,      # R in angstrem
                    x,      # xd
                    x * m,  # xm
                    1,      # eps_d
                    m * m,  # eps_m
                    d_parl[i],      # d_parallel
                    d_perp[i])      # d_perp
    mesomie.calc_Q()
    Qext.append(mesomie.GetQext())

    mie.SetLayersSize(x)
    mie.SetLayersIndex(m)
    mie.RunMieCalculation()
    Qext_mie.append(mie.GetQext())

    # m = m_drude
    # mesomie.calc_ab(R*10,      # R in angstrem
    #                 x,      # xd
    #                 x * m,  # xm
    #                 1,      # eps_d
    #                 m * m,  # eps_m
    #                 d_parl[i],      # d_parallel
    #                 d_perp[i])      # d_perp
    # mesomie.calc_Q()
    # Qext_drude_nc.append(mesomie.GetQext())

    # print(x, m, Qext[-1] - mie.GetQext())

    om_rat_plot.append(om_star_rat)

plt.plot(om_rat_plot, Qext_mie,
         label='classic', color='black', lw=4)

# plt.plot(om_rat_plot, np.real(eps_inf_drude),
#          label='real drude', color='blue', lw=1)
# plt.plot(om_rat_plot, np.imag(eps_inf_drude),
#          label='imag drude', color='red', lw=1)

# plt.plot(om_rat_plot, np.real(eps_m_drude),
#          label='real drude', color='blue', lw=2)
# plt.plot(om_rat_plot, np.imag(eps_m_drude),
#          label='imag drude', color='red', lw=2)


plt.plot(om_rat_plot, Qext,
         label='non-classic', color='red', lw=4)
# plt.plot(om_rat_plot, Qext_drude_nc,
#          label='non-classic drude fixed\nomega_p = 9.02 eV\ngamma = 0.15eV\neps_inf_drud = 4.65', color='blue', lw=2)
for j in range(7):
    Qext_drude = []
    step = 0.02
    for i in range(len(omega_star_ratio)):
        om_star_rat = omega_star_ratio[i]
        if (om_star_rat < min_lim_omega_star_ratio
                or om_star_rat > max_lim_omega_star_ratio):
            continue
        omega = om_star_rat*omega_p_star
        WL_mkm = 2*np.pi/((omega/c)/h_reduced)*1e6
        if WL_mkm < min_WL_available or WL_mkm > max_WL_available:
            continue
        x_const = (omega/c) * 1e-9/h_reduced

        x_cs = [(R-shell_h) * x_const, R * x_const],
        m = index_Ag[i, 1]
        m_drude = cmath.sqrt(eps_m(omega, eps_inf(omega, m**2), omega_p))
        m_drude_shell = cmath.sqrt(
            eps_m(omega, eps_inf(omega, m**2), omega_p*(0.96+step*j)))
        m_cs = [m_drude, m_drude_shell]
        mie.SetLayersSize(x_cs)
        mie.SetLayersIndex(m_cs)
        mie.RunMieCalculation()
        Qext_drude.append(mie.GetQext())
    plt.plot(om_rat_plot, Qext_drude,
             label=f'omega_p*{((0.96+step*j)*100)/100}', color='gray', lw=2)

plt.legend()
# plt.yscale('log')
plt.xlim((min_lim_omega_star_ratio, max_lim_omega_star_ratio))
# plt.ylim((y_min, y_max))
plt.title(
    "R="+str(R)+f'\nfor core-shell totalR is the same,\nshell_h={shell_h}')
plt.show()
