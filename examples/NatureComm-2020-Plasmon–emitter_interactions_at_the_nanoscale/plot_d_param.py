#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
from cProfile import label
from matplotlib import markers, pyplot as plt
import numpy as np
from scipy import interpolate

arr2D = np.loadtxt('rs4-im-d_perp.csv', delimiter=',')
im_d = arr2D[arr2D[:, 0].argsort()]
arr2D = np.loadtxt('rs4-re-d_perp.csv', delimiter=',')
re_d = arr2D[arr2D[:, 0].argsort()]

xmin_im = np.min(im_d[:, 0])
xmin_re = np.min(re_d[:, 0])
xmax_im = np.max(im_d[:, 0])
xmax_re = np.max(re_d[:, 0])
x = np.linspace(np.max([xmin_im, xmin_re]), np.min([xmax_im, xmax_re]), 1000)
im_d_y = interpolate.interp1d(im_d[:, 0],  im_d[:, 1])
re_d_y = interpolate.interp1d(re_d[:, 0],  re_d[:, 1])

data = np.array([x.T, re_d_y(x).T, im_d_y(x).T])
np.savetxt('rs4-d_perp.txt', data)

# print(data)
# plt.plot(im_d[:, 0], im_d[:, 1], marker='o', ls='')
# plt.plot(x, im_d_y(x))
# plt.plot(re_d[:, 0], re_d[:, 1], marker='o', ls='')
# plt.plot(x, re_d_y(x))
# plt.xlim((0, 1))
# plt.ylim((-4, 5))
# plt.show()

from_disk = np.loadtxt('rs4-d_perp.txt')
plt.plot(from_disk[0, :], from_disk[1, :], label='re d')
plt.plot(from_disk[0, :], from_disk[2, :], label='im d')
plt.legend()
plt.show()
