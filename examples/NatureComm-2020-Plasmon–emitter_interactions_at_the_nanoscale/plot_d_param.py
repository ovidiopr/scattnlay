#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
from cProfile import label
from matplotlib import markers, pyplot as plt
import numpy as np
from scipy import interpolate

import scipy.io
mat = scipy.io.loadmat('d-parameters/rs=4.mat')
x_mat = mat['omegav'][0]
d_perp_mat = mat['dperp'][0]*10
x = np.linspace(x_mat[0], x_mat[-1], 1001)
im_d_y = interpolate.interp1d(x_mat,  np.imag(d_perp_mat), kind='cubic')
re_d_y = interpolate.interp1d(x_mat,  np.real(d_perp_mat))
data = np.array([x.T, re_d_y(x).T, im_d_y(x).T])
np.savetxt('rs4-d_perp_interpolated.txt', data)

mat = scipy.io.loadmat('d-parameters/silver_xd_1.18A.mat')
x_mat = mat['omegav'][0]
d_perp_mat = mat['dperp'][0]*10
d_parl_mat = mat['dparl'][0]*10
x = np.linspace(x_mat[0], x_mat[-1], 1001)
kind = 'linear'
im_d_y = interpolate.interp1d(x_mat,  np.imag(d_perp_mat), kind=kind)
re_d_y = interpolate.interp1d(x_mat,  np.real(d_perp_mat))
data = np.array([x.T, re_d_y(x).T, im_d_y(x).T])
np.savetxt('silver-d_perp_interpolated.txt', data)
im_d_y = interpolate.interp1d(x_mat,  np.imag(d_parl_mat), kind=kind)
re_d_y = interpolate.interp1d(x_mat,  np.real(d_parl_mat))
data = np.array([x.T, re_d_y(x).T, im_d_y(x).T])
np.savetxt('silver-d_parl_interpolated.txt', data)


from_disk = np.loadtxt('rs4-d_perp_interpolated.txt')
plt.figure('rs4')
plt.plot(from_disk[0, :], from_disk[1, :], label='re d')
plt.plot(from_disk[0, :], from_disk[2, :], label='im d')
plt.plot(x_mat, np.real(d_perp_mat), label='re d mat')
plt.plot(x_mat, np.imag(d_perp_mat), label='re d mat')
plt.legend()
plt.show()


from_disk = np.loadtxt('silver-d_perp_interpolated.txt')
plt.figure('silver')
plt.plot(from_disk[0, :], from_disk[1, :], label='re d perp')
plt.plot(from_disk[0, :], from_disk[2, :], label='im d perp')
from_disk = np.loadtxt('silver-d_parl_interpolated.txt')
plt.plot(from_disk[0, :], from_disk[1, :], label='re d parl')
plt.plot(from_disk[0, :], from_disk[2, :], label='im d parl')

plt.legend()
plt.show()
