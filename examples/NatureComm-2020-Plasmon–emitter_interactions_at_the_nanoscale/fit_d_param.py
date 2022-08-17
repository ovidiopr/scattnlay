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
d_perp_rs4 = from_disk[1, ::step] + 1j*from_disk[2, ::step]


def rms(x0):
    d_fit = spectra(omega_ratio, x0)
    diff_re = np.real(d_perp_rs4 - d_fit)
    rms = np.sqrt(np.sum(np.abs(diff_re)**2))
    diff_im = np.imag(d_perp_rs4 - d_fit)
    rms += np.sqrt(np.sum(np.abs(diff_im)**2))
    return rms


poles = 1
dim = poles*4
x0 = np.random.random(dim)


problem_dict1 = {
    "fit_func": rms,
    "lb": [0,  0, -10, -10],
    "ub": [10, 10, 10, 10],
    "minmax": "min",
}
epoch = 300
pop_size = 10
model = AdaptiveEO(problem_dict1, epoch, pop_size)
# x, best_fitness = model.solve()
# print(f"Solution: {x}, Fitness: {best_fitness}")

x, es = cma.fmin2(rms, x0, sigma0=0.2)


d_fit = spectra(omega_ratio, x)
print('first round x =', x)


print('x =', x)
# print('ex =', es)

# print('d_fit =', d_fit)

plt.figure('rs4')
# plt.title('rms = '+str(rms(x)/x.size))
plt.plot(omega_ratio, np.real(d_perp_rs4), label='re d')
plt.plot(omega_ratio, np.imag(d_perp_rs4), label='im d')

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
