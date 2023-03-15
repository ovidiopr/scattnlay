#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
from glob import glob
from operator import truediv
from matplotlib import pyplot as plt
import cmath
from scattnlay import mesomie, mie
import numpy as np
import scipy.io
from functools import wraps
import time
from eval_spectra import multi_lorentzian, params_count as pc
import cma
import sys
from optical_constants import read_refractive_index_from_yaml as read_nk

isPlot = True
# isPlot = False


def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(
            f'Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds')
        return result
    return timeit_wrapper


min_lim = 0.4
max_lim = 0.75

# mat = scipy.io.loadmat('d-parameters/rs=4.mat')
# omega_ratio = mat['omegav'][0]
# d_perp = mat['dperp'][0]*10

from_disk = np.loadtxt('rs4-d_perp_interpolated.txt')
from_disk = from_disk[:,
                      from_disk[0, :] > min_lim
                      ]
from_disk = from_disk[:,
                      from_disk[0, :] < max_lim
                      ]

omega_ratio = from_disk[0, :]
d_perp = from_disk[1, :] + 1j*from_disk[2, :]

c = 299792458  # m/s
h_reduced = 6.5821e-16  # eV s
omega_p = 5.9  # eV
gamma = 0.1  # eV
eps_d = 1

# omega_p_star = 3.81  # eV #orig
omega_p_star = 5.81  # eV
omega = omega_ratio*omega_p_star
# WL_mkm = 2*np.pi/((omega/c)/h_reduced)*1e6

WL = 2*np.pi/((omega_ratio * omega_p_star/c)/h_reduced)*1e6  # mkm
# min_WL_available = 0.1879
# max_WL_available = 1.9370
# WL[WL < min_WL_available] = min_WL_available
# WL[WL > max_WL_available] = max_WL_available
index_Ag = read_nk('Ag-Johnson-1972.yml', WL, kind=1)
eps_Ag = index_Ag**2
# print(index_Ag)
# plt.figure('index')
# plt.plot(index_Ag[:, 0], np.real(index_Ag[:, 1]))
# plt.plot(index_Ag[:, 0], np.imag(index_Ag[:, 1]))
# plt.show()
# sys.exit()


def eps_m(omega):
    return 1 - omega_p * omega_p / (omega*omega + 1j*omega*gamma)


n = 1 if isPlot else 1
Rs = [2.5]*n  # nm
y_min = [1e-2]*n
y_max = [1e1]*n

# Rs = [2.5, 5, 10, 25]  # nm
# y_min = [1e-2, 1e-2, 1e-1, 1e-1]
# y_max = [1e1, 1e1, 5e1, 5e1]

# for om_rat in omega_ratio:
R = 0


@ timeit
def run():
    for fig in range(len(Rs)):
        global R
        R = Rs[fig]
        Qext = []
        Qext_fit = []
        om_rat_plot = []
        x_fit_from_d_perp = np.array([0.09391149,  0.81234806, -0.30526326, -0.01044856,
                                      0.3121473,   0.55333457, -0.01684191,  0.04727765,
                                      0.11249052,  0.88368699, -0.04680872, -0.10982548])
        omega = omega_ratio*omega_p
        x = (omega/c) * R * 1e-9/h_reduced
        m_orig = np.sqrt(eps_m(omega))
        # m = np.sqrt(eps_m(omega))

        # m = index_Ag
        m = index_Ag[:, 1]
        # plt.figure('index2')
        # plt.plot(x, np.real(m), label='orig')
        # plt.plot(x, np.imag(m), label='orig')
        # plt.plot(x, np.real(m), label='mod')
        # plt.plot(x, np.imag(m), label='mod')
        # plt.legend()
        # plt.show()
        # sys.exit()

        Qext = getMesoQext(R, x, m_orig, d_perp)

        def getQfit(x_fit):
            d_fit = multi_lorentzian(omega_ratio, x_fit)
            return getMesoQext(R, x, m, d_fit)

        def rms_base(Q_rms, Q_fit):
            diff_re = np.real(Q_rms - Q_fit)
            diff_im = np.imag(Q_rms - Q_fit)
            rms = (np.sum(np.abs(diff_re)**2))
            rms += (np.sum(np.abs(diff_im)**2))
            return rms

        def rms(x0):
            Q_fit = getQfit(x0)
            return rms_base(Q_rms, Q_fit)

        def rms_x12(x0):
            x0[:pc] = x0[:pc]+x1
            Q_fit = getQfit(x0)
            return rms_base(Q_rms, Q_fit)

        def rms_x1(x0):
            Q_fit = getQfit(np.hstack((x1, x0)))
            return rms_base(Q_rms, Q_fit)

        x_fit_from_d_perp = np.array([0.09391149,  0.81234806, -0.30526326, -0.01044856,
                                      0.3121473,   0.55333457, -0.01684191,  0.04727765,
                                      0.11249052,  0.88368699, -0.04680872, -0.10982548])

        x0 = np.random.random(pc)
        Q_rms = Qext

        xfit, es = cma.fmin2(rms, x0, sigma0=2)
        # xfit = np.array([0.32740616895990493, -0.844352181880013, -
        #                  0.5682552466755, 0.015925488550861087])
        x1 = np.copy(xfit)
        Qext_fit = getQfit(xfit)
        # print(x1)

        # x0 = np.random.random(pc)
        # # x0 = np.copy(x1)
        # x0 = np.hstack((x1, x0))
        # xfit, es = cma.fmin2(rms, x0, sigma0=0.2)
        # # xfit = np.array([-9.03308419e-01, -1.77584180e+00, -2.36259518e+01, 3.41837008e+00,
        #                  5.75724471e-06, , 2.64230613e+00, , 4.72267415e+01, -1.61624064e+00])
        # # xfit = np.array([0.32740616895990493, -0.844352181880013, -
        # #  0.5682552466755, 0.015925488550861087])
        # x2 = np.copy(xfit)
        # Qext_fit = getQfit(np.hstack((x1, x2)))

        print(xfit)

        om_rat_plot = omega_ratio

        if isPlot:
            Qext_fit = getQfit(xfit)
            Qext_no_d = getMesoQext(R, x, m_orig, np.zeros(len(x)))
            Qext_no_d_mod = getMesoQext(R, x, m, np.zeros(len(x)))
            d_fit = multi_lorentzian(omega_ratio, xfit)
            plt.figure(fig)
            plt.plot(om_rat_plot, Qext_fit,
                     label='fitted orig using mod', color='gray', lw=4)
            plt.plot(om_rat_plot, Qext_no_d,
                     label='d=0 orig', color='gray', lw=4)
            plt.plot(om_rat_plot, Qext_no_d_mod,
                     label='d=0 mod', color='gray', lw=4)
            plt.plot(om_rat_plot, Qext, label='direct orig', color='red', lw=1)
            plt.legend()
            # plt.yscale('log')
            # plt.xlim((0.4, 0.404))
            # plt.ylim((y_min[fig]*1.9, y_min[fig]*2.3))
            plt.xlim((0.4, 0.75))
            # plt.xlim((0.5, 0.6))
            # plt.ylim((y_min[fig], y_max[fig]))
            plt.title("Qext for R="+str(R))

            plt.figure("d_param "+str(fig))
            plt.plot(om_rat_plot, np.real(d_fit), label='re fit', lw=4)
            plt.plot(om_rat_plot, np.real(d_perp), label='re d_perp')
            plt.plot(om_rat_plot, np.imag(d_fit), label='im fit', lw=4)
            plt.plot(om_rat_plot, np.imag(d_perp), label='im d_perp')
            plt.legend()

    if isPlot:
        plt.show()


# @timeit
def getMesoQext(R, x, m, d_perp):
    Qext = np.zeros((len(x)))
    for i in range(len(x)):
        mesomie.calc_ab(R*10,      # calc_ab needs R in angstrem
                        x[i],      # xd
                        x[i] * m[i],  # xm
                        1,      # eps_d
                        m[i] * m[i],  # eps_m
                        0,      # d_parallel
                        d_perp[i])      # d_perp
        mesomie.calc_Q()
        Qext[i] = mesomie.GetQext()
    return Qext


run()
