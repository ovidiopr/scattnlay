#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
from matplotlib import pyplot as plt
import cmath
from evalMie import *
import numpy as np
min_lim = 0.4
max_lim = 0.75

# mat = scipy.io.loadmat('d-parameters/rs=4.mat')
# omega_ratio = mat['omegav'][0]
# d_perp = mat['dperp'][0]*10

from_disk = np.loadtxt('rs4-d_perp_interpolated.txt')
omega_ratio = from_disk[0, :]
d_perp = from_disk[1, :] + 1j*from_disk[2, :]

c = 299792458  # m/s
h_reduced = 6.5821e-16  # eV s
omega_p = 5.9  # eV
gamma = 0.1  # eV
eps_d = 1


def eps_m(omega):
    return 1 - omega_p * omega_p / (omega*omega + 1j*omega*gamma)




Rs = [2.5, 5, 10, 25]  # nm
y_min = [1e-2, 1e-2, 1e-1, 1e-1]
y_max = [1e1, 5e1, 5e1, 5e1]

# for om_rat in omega_ratio:
for fig in range(len(Rs)):
    R = Rs[fig]
    Qext = []
    Qext_mie = []
    om_rat_plot = []
    x_in = []
    m_in = []
    d_perp_in = []
    for i in range(len(omega_ratio)):
        om_rat = omega_ratio[i]
        if om_rat < min_lim or om_rat > max_lim:
            continue
        omega = om_rat*omega_p
        m = cmath.sqrt(eps_m(omega))
        x = (omega/c) * R * 1e-9/h_reduced
        x_in.append(x)
        m_in.append(m)
        om_rat_plot.append(om_rat)
        d_perp_in.append(d_perp[i])

    Qext_mie, _ = eval_mie_spectrum(x_in, m_in)
    Qext, _ = eval_mesomie_spectrum(x_in, m_in, R, d_perp_in)

    print(Qext)
    plt.figure(fig)
    plt.plot(om_rat_plot, Qext_mie, label='classic', color='gray', lw=4)
    plt.plot(om_rat_plot, Qext, label='non-classic', color='red', lw=4)
    # plt.plot(om_rat_plot, Qext, label="R="+str(R), lw=1)
    plt.legend()
    plt.yscale('log')
    plt.xlim((0.4, 0.75))
    plt.ylim((y_min[fig], y_max[fig]))
    plt.title("R="+str(R))
plt.show()
