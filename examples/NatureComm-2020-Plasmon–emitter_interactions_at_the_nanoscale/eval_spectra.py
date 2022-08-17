from numba import njit, prange, complex128, float64
import numpy as np


@njit(complex128(float64, float64[:]),
      fastmath=True, cache=True)
def lorentzian(omega, xvec):
    res = 0
    poles = len(xvec)//4
    for i in range(poles):
        gamma = xvec[i+0]
        omega_n = xvec[i+1]
        f = xvec[i+2] + 1j*xvec[i+3]
        res = res + f / (omega * (omega + 1j*gamma) - omega_n**2)
    return res


@njit(complex128[:](float64[:], float64[:]),
      parallel=True, fastmath=True, cache=True)
def spectra(omega, xvec):
    poles = len(xvec)//4
    for i in range(poles):
        xvec[i] = np.absolute(xvec[i])
        xvec[i+1] = np.absolute(xvec[i+1])

    val = np.zeros(omega.size, dtype=np.cdouble)
    for i in prange(omega.size):
        val[i] = lorentzian(omega[i], xvec)
    return val
