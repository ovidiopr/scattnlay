from ast import Constant
from numba import njit, prange, complex128, float64
import numpy as np


# omega_init = True
params_count = 4


@njit(complex128(float64, float64[:]),
      fastmath=True, cache=True)
def lorentzian(omega, xvec):
    # global omega_init
    pc = params_count
    res = 0
    poles = len(xvec)//params_count
    factor = 1
    for i in range(poles):
        gamma = xvec[pc*i+0] / (factor**i)
        omega_n = xvec[pc*i+1] / (factor**i)
        f = (xvec[pc*i+2] + 1j*xvec[pc*i+3]) / (factor**i)
        # if (omega_init):
        #     print('init', gamma, omega_n, f)
        if (np.abs(f) == 0):
            return res
        res = res + f / (omega * (omega + 1j*gamma) - omega_n**2)
    # omega_init = False
    return res


# @njit(complex128[:](float64[:], float64[:]),
#       parallel=True, fastmath=True, cache=True)
def spectra(omega, xvec):
    # global omega_init
    # omega_init = True
    xvec_in = xvec
    poles = len(xvec)//params_count
    # for i in range(poles):
    #     for j in range(2):
    #         if (xvec_in[i+j] < 0):
    #             xvec_in[i+j] = 0
    #     # if (xvec_in[i+2] > 0):
    #         # xvec_in[i+2] = 0

    val = np.zeros(omega.size, dtype=np.cdouble)
    # print('in eval', xvec_in)
    for i in range(omega.size):
        val[i] = lorentzian(omega[i], xvec_in)
    return val
