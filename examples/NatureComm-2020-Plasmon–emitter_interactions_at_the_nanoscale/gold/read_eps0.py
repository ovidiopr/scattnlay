# -*- coding: utf-8 -*-
# GPL license from Smuthi project

from scipy.interpolate import interp1d
import numpy as np
from matplotlib import pyplot as plt


def read_eps0(filename, vacuum_wavelength, kind=1):
    data_np = np.loadtxt(filename, comments='#', delimiter=',')
    data_wl = data_np[:, 0]
    eps0 = data_np[:, 1]
    f = interp1d(data_wl, eps0, kind=kind)
    data_out = np.transpose(np.vstack(
        (vacuum_wavelength, f(vacuum_wavelength))))
    if len(data_out) == 1:
        return data_out[0]
    return data_out


WLs = np.linspace(450, 650, 1001)
eps0_real = read_eps0('Au_eps0_real.csv', WLs, kind=2)
eps0_real1 = read_eps0('Au_eps0_real.csv', WLs, kind=1)
eps0_imag = read_eps0('Au_eps0_imag.csv', WLs)
plt.plot(eps0_real[:, 0], eps0_real[:, 1])
plt.plot(eps0_real1[:, 0], eps0_real1[:, 1])
plt.plot(eps0_imag[:, 0], eps0_imag[:, 1])
plt.show()
