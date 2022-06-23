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

# arr2D = np.loadtxt('rs4-im-d_perp.csv', delimiter=',')
# im_d = arr2D[arr2D[:, 0].argsort()]
# arr2D = np.loadtxt('rs4-re-d_perp.csv', delimiter=',')
# re_d = arr2D[arr2D[:, 0].argsort()]

x = np.linspace(x_mat[0], x_mat[-1], 1001)
im_d_y = interpolate.interp1d(x_mat,  np.imag(d_perp_mat), kind='cubic')
re_d_y = interpolate.interp1d(x_mat,  np.real(d_perp_mat))

data = np.array([x.T, re_d_y(x).T, im_d_y(x).T])
np.savetxt('rs4-d_perp_interpolated.txt', data)


from_disk = np.loadtxt('rs4-d_perp_interpolated.txt')
plt.plot(from_disk[0, :], from_disk[1, :], label='re d')
plt.plot(from_disk[0, :], from_disk[2, :], label='im d')
plt.plot(x_mat, np.real(d_perp_mat), label='re d mat')
plt.plot(x_mat, np.imag(d_perp_mat), label='re d mat')
plt.legend()
plt.show()
