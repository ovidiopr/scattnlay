from ast import Constant
from numba import njit, prange, complex128, float64
import numpy as np
from scattnlay import mesomie

params_count = 4


@njit(complex128(float64, float64[:]),
      fastmath=True, cache=True)
def lorentzian(omega, xvec):
    pc = params_count
    res = 0
    poles = len(xvec)//params_count
    for i in range(poles):
        gamma = np.abs(xvec[pc*i+0])  # >0
        omega_n = np.abs(xvec[pc*i+1])  # >0
        f = (xvec[pc*i+2] + 1j*xvec[pc*i+3])
        if (np.abs(f) == 0):
            return res
        res = res + f / (omega * (omega + 1j*gamma) - omega_n**2)
    return res


@njit(complex128[:](float64[:], float64[:]),
      fastmath=True, cache=True)
def multi_lorentzian(omega, xvec):
    poles = len(xvec)//params_count
    val = np.zeros(omega.size, dtype=np.cdouble)
    for i in range(omega.size):
        val[i] = lorentzian(omega[i], xvec)
    return val
