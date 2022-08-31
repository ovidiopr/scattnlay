#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
from matplotlib import pyplot as plt
import numpy as np
import cma
from eval_spectra import multi_lorentzian, params_count as pc
import sys

from_disk = np.loadtxt('rs4-d_perp_interpolated.txt')
# from_disk = from_disk[:,
#                       from_disk[0, :] > 0.4
#                       ]
# from_disk = from_disk[:,
#                       from_disk[0, :] < 0.8
#                       ]

step = 5
omega_ratio = np.copy(from_disk[0, ::step])
d_rs4 = from_disk[1, ::step] + 1j*from_disk[2, ::step]


def rms(x0):
    d_fit = multi_lorentzian(omega_ratio, x0)
    diff_re = np.real(d_rms - d_fit)
    diff_im = np.imag(d_rms - d_fit)
    # rms = np.sqrt(np.sum(np.abs(diff_re)**2))
    # rms += np.sqrt(np.sum(np.abs(diff_im)**2))
    rms = (np.sum(np.abs(diff_re)**2))
    rms += (np.sum(np.abs(diff_im)**2))
    return rms


x0 = np.random.random(pc)
d_rms = d_rs4
x, es = cma.fmin2(rms, x0, sigma0=0.2)
# x = np.array([0.1332754793043937, 0.8248539955310855, -
#              0.5024906405674647, -0.08076797734842008])
x1 = np.copy(x)


x0 = np.random.random(pc)
d_rms = d_rs4 - multi_lorentzian(omega_ratio, x1)
x, es = cma.fmin2(rms, x0, sigma0=2)
# x = np.array([0.0019369902204593222, 0.9978752165162739,
#              0.05801769075873917, -0.032110084386726336])
x2 = np.copy(x)


x0 = np.hstack((x1, x2))
d_rms = d_rs4
x, es = cma.fmin2(rms, x0, sigma0=0.02)
# x = np.array([0.11878109939932953, 0.8142072049467617, -0.43466301805510305, -0.012772472816983446,
#               0.012573279034397847, 1.0010563127596699, 0.07665592968844606, -0.07679477702750433])
x12 = x


x0 = np.random.random(pc)
d_rms = d_rs4 - multi_lorentzian(omega_ratio, x12)
x, es = cma.fmin2(rms, x0, sigma0=0.2)
# x = np.array([0.1434266344187499, 0.5188802822956783, -
#              0.00950433613285183, 0.013585684987833985])
x3 = np.copy(x)


x0 = np.hstack((x1, x2, x3))
d_rms = d_rs4
x, es = cma.fmin2(rms, x0, sigma0=0.02)
# x = np.array([0.09914803, 0.81158511, -0.34941712, 0.01388308,
#               0.01560184, 1.00551237, 0.11006553, -0.09818891,
#               0.34432793, 0.64182428, -0.0803532, 0.04641341])
x123 = x

#               [ 0.09349663 -0.8188804  -0.36874431 -0.08079194
#                 0.25712096  0.56012451 -0.02089828  0.03133694]

#               [ 0.09391149  0.81234806 -0.30526326 -0.01044856
#                 0.3121473   0.55333457 -0.01684191  0.04727765
#                 0.11249052  0.88368699 -0.04680872 -0.10982548]
print(x)

d_fit = multi_lorentzian(omega_ratio, x)
plt.figure('rs4')
plt.plot(omega_ratio, np.real(d_rms), label='re d')
plt.plot(omega_ratio, np.imag(d_rms), label='im d')
plt.plot(omega_ratio, np.real(d_fit), label='re d fit', alpha=0.2, lw=3)
plt.plot(omega_ratio, np.imag(d_fit), label='im d fit', alpha=0.2, lw=3)
plt.axhline(y=0.0, color='black', linestyle='-', lw=1, alpha=0.2)
plt.legend()
plt.show()
