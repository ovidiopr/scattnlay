#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
from functools import lru_cache
from matplotlib import markers, pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import cma
# import pyfde
from numba import njit, float64
from eval_spectra import spectra

from mealpy.physics_based.EO import AdaptiveEO


from_disk = np.loadtxt('rs4-d_perp_interpolated.txt')
step = 5
omega_ratio = np.copy(from_disk[0, ::step])
d_rs4 = from_disk[1, ::step] + 1j*from_disk[2, ::step]

d_rms = d_rs4


def rms(x0):
    d_fit = spectra(omega_ratio, x0)
    diff_re = np.real(d_rms - d_fit)
    rms = np.sqrt(np.sum(np.abs(diff_re)**2))
    diff_im = np.imag(d_rms - d_fit)
    rms += np.sqrt(np.sum(np.abs(diff_im)**2))
    return rms


poles = 1
dim = poles*4
x0 = np.random.random(dim)
d_rms = d_rs4
# x, es = cma.fmin2(rms, x0, sigma0=0.2)
x = np.array([0.13421489, 0.82250415, -0.50359304, -0.0591722])
x1 = x

x0 = np.random.random(dim)
d_rms = d_rs4 - spectra(omega_ratio, x1)
x, es = cma.fmin2(rms, x0, sigma0=2)
# x = [0.00051888486267353, 0.9991499897780783,
#      0.06926055304806265, -0.030608812209114836] # fitness = 4.7
x2 = x

d_fit = spectra(omega_ratio, x)

plt.figure('rs4')
# plt.title('rms = '+str(rms(x)/x.size))
plt.plot(omega_ratio, np.real(d_rms), label='re d')
plt.plot(omega_ratio, np.imag(d_rms), label='im d')

# plt.plot(omega_ratio, np.real(
#     d_rs4 - spectra(omega_ratio, x1)), label='diff re d')
# plt.plot(omega_ratio, np.imag(
#     d_rs4 - spectra(omega_ratio, x1)), label='diff im d')

plt.plot(omega_ratio, np.real(d_fit), label='re d fit', alpha=0.2, lw=3)
plt.plot(omega_ratio, np.imag(d_fit), label='im d fit', alpha=0.2, lw=3)

# plt.plot(omega_ratio, func(xdata, *popt), 'r-',

#          label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
plt.legend()
plt.show()


# from_disk = np.loadtxt('silver-d_perp_interpolated.txt')
# plt.figure('silver')
# plt.plot(from_disk[0, :], from_disk[1, :], label='re d perp')
# plt.plot(from_disk[0, :], from_disk[2, :], label='im d perp')
# from_disk = np.loadtxt('silver-d_parl_interpolated.txt')
# plt.plot(from_disk[0, :], from_disk[1, :], label='re d parl')
# plt.plot(from_disk[0, :], from_disk[2, :], label='im d parl')

# plt.legend()
# plt.show()
