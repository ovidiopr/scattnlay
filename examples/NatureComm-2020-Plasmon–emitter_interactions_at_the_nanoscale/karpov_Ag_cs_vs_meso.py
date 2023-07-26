#!/usr/bin/env python3
import cmath
import numpy as np
from matplotlib import pyplot as plt

from optical_constants import read_refractive_index_from_yaml as read_nk

from evalMie import *
from karpov_materials import *

WL_min = 325 #nm
WL_max = 390 #nm

# WL_min = 335 #nm
# WL_max = 500 #nm
shell_h = 0.5
n_m = 1

d_param_filename = ['silver-d_perp_interpolated.txt',
                    'silver-d_parl_interpolated.txt']


def eval_xm(WLs,R, shell_h, n_m):
    x_cs_in, m_cs_in = [], []
    x_bulk_in, m_bulk_in = [], []

    w,wp, Gamma = evalDrudeParams(R, WLs)
    
    for i in range(len(WLs)):
        x_const = 2 * n_m *np.pi/ WLs[i]
        x_cs = [(R-shell_h) * x_const, R * x_const],
        x_meso = R*x_const 

        n_drude = cmath.sqrt(
            eps_drude(w[i], wp, Gamma))
        n_shell = cmath.sqrt(
            eps_drude(w[i], wp * correction_wp[j], Gamma))
        m_cs = [n_drude / n_m, n_shell / n_m]
        m_meso = n_drude / n_m

        x_cs_in.append(x_cs)
        m_cs_in.append(m_cs)
        x_bulk_in.append(x_meso)
        m_bulk_in.append(m_meso)
    return x_cs_in, m_cs_in, x_bulk_in, m_bulk_in

def plot_Q(WLs, Qext_drude, legend, ls='-'):
    plt.plot(WLs, Qext_drude, ls=ls)
    plt.legend(legend, fontsize=14)
    plt.xlabel("\u03BB, nm", fontsize=14)
    plt.ylabel("Normalized extintion cross-section, $nm^{2}$", fontsize=14)


# correction_wp = [0.97, 0.95, 0.92, 0.89, 0.87, 0.855, 0.84, 0.83]

correction_wp = [0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.9]
Rs = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]

correction_wp = correction_wp[:3]
Rs = Rs[:3]


legend_cs = [str(r)+'-'+str(w)+' cs' for r,w in zip(Rs, correction_wp)]
legend_mie = [str(r)+'-'+str(w) for r,w in zip(Rs, correction_wp)]
legend_meso = [str(r)+'-'+str(w)+' d' for r,w in zip(Rs, correction_wp)]

WLs_nm, d_perp_in, d_parl_in = get_d_params(WL_min, WL_max, d_param_filename)

index_Ag = read_nk('Ag-Johnson-1972.yml', WLs_nm, units='nm', kind=1)
m_Johnson = index_Ag[:, 1]

# WLs_nm = np.linspace(WL_min, WL_max, 401)
Qext = []
for j in range(len(Rs)):

    x_cs_in, m_cs_in, x_bulk_in, m_bulk_in = eval_xm(WLs_nm,Rs[j], shell_h, n_m)

    Qext_cs_drude, _ = eval_mie_spectrum(x_cs_in, m_cs_in)
    plt.figure('mie')
    plot_Q(WLs_nm, Qext_cs_drude, legend_cs, ls='dotted')

    # m_bulk_in = m_Johnson

    Qext_bulk_drude, _ = eval_mie_spectrum(x_bulk_in, m_bulk_in)
    plt.figure('mie')
    plot_Q(WLs_nm, Qext_bulk_drude, legend_mie)

    Qext_mesomie, _ = eval_mesomie_spectrum(x_bulk_in, m_bulk_in, Rs[j], d_perp_in, d_parl_in)
    plt.figure('mie')
    plot_Q(WLs_nm, Qext_mesomie, legend_meso, ls='--')

    # print(WLs_nm, Qext_drude)
    # print(max(Qext_drude))

plt.show()

