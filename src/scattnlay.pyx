#    Copyright (C) 2009-2018 Ovidio Pena <ovidio@bytesfall.com>
#    Copyright (C) 2013-2018 Konstantin Ladutenko <kostyfisik@gmail.com>
#
#    This file is part of scattnlay
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    The only additional remark is that we expect that all publications
#    describing work using this software, or all commercial products
#    using it, cite at least one of the following references:
#    [1] O. Peña and U. Pal, "Scattering of electromagnetic radiation by
#        a multilayered sphere," Computer Physics Communications,
#        vol. 180, Nov. 2009, pp. 2348-2354.
#    [2] K. Ladutenko, U. Pal, A. Rivera, and O. Peña-Rodríguez, "Mie
#        calculation of electromagnetic near-field for a multilayered
#        sphere," Computer Physics Communications, vol. 214, May 2017,
#        pp. 225-230.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# distutils: language = c++
# distutils: sources = nmie.cc

from __future__ import division
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.vector cimport complex

cdef inline double *npy2c(np.ndarray a):
    assert a.dtype == np.float64

    if not (<object>a).flags["C_CONTIGUOUS"]: # Array is not contiguous, need to make contiguous copy
        a = a.copy('C')

    # Return data pointer
    return <double *>(a.data)

cdef extern from "py_nmie.h":
    cdef int ScattCoeffs(int L, int pl, vector[double] x, vector[complex] m, int nmax, double anr[], double ani[], double bnr[], double bni[])
    cdef int nMie(int L, int pl, vector[double] x, vector[complex] m, int nTheta, vector[double] Theta, int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, double S1r[], double S1i[], double S2r[], double S2i[])
    cdef int nField(int L, int pl, vector[double] x, vector[complex] m, int nmax, int nCoords, vector[double] Xp, vector[double] Yp, vector[double] Zp, double Erx[], double Ery[], double Erz[], double Eix[], double Eiy[], double Eiz[], double Hrx[], double Hry[], double Hrz[], double Hix[], double Hiy[], double Hiz[])

def scattcoeffs(np.ndarray[np.float64_t, ndim = 2] x, np.ndarray[np.complex128_t, ndim = 2] m, np.int_t nmax, np.int_t pl = -1):
    """
    scattcoeffs(x, m[, nmax, pl])

    Calculate the scattering coefficients required to calculate both the
    near- and far-field parameters.

        x: Size parameters (2D ndarray)
        m: Relative refractive indices (2D ndarray)
        nmax: Maximum number of multipolar expansion terms to be used for the
              calculations. Only use it if you know what you are doing, otherwise
              set this parameter to -1 and the function will calculate it.
        pl: Index of PEC layer. If there is none just send -1

    Returns: (terms, an, bn)
    with
        terms: Number of multipolar expansion terms used for the calculations
        an, bn: Complex scattering coefficients
    """
    cdef Py_ssize_t i

    cdef np.ndarray[np.int_t, ndim = 1] terms = np.zeros(x.shape[0], dtype = np.int)

    cdef np.ndarray[np.complex128_t, ndim = 2] an = np.zeros((x.shape[0], nmax), dtype = np.complex128)
    cdef np.ndarray[np.complex128_t, ndim = 2] bn = np.zeros((x.shape[0], nmax), dtype = np.complex128)

    cdef np.ndarray[np.float64_t, ndim = 1] anr
    cdef np.ndarray[np.float64_t, ndim = 1] ani
    cdef np.ndarray[np.float64_t, ndim = 1] bnr
    cdef np.ndarray[np.float64_t, ndim = 1] bni

    for i in range(x.shape[0]):
        anr = np.zeros(nmax, dtype = np.float64)
        ani = np.zeros(nmax, dtype = np.float64)
        bnr = np.zeros(nmax, dtype = np.float64)
        bni = np.zeros(nmax, dtype = np.float64)

        terms[i] = ScattCoeffs(x.shape[1], pl, x[i].copy('C'), m[i].copy('C'), nmax, npy2c(anr), npy2c(ani), npy2c(bnr), npy2c(bni))

        an[i] = anr.copy('C') + 1.0j*ani.copy('C')
        bn[i] = bnr.copy('C') + 1.0j*bni.copy('C')

    return terms, an, bn

def scattnlay(np.ndarray[np.float64_t, ndim = 2] x, np.ndarray[np.complex128_t, ndim = 2] m, np.ndarray[np.float64_t, ndim = 1] theta = np.zeros(0, dtype = np.float64), np.int_t nmax = -1, np.int_t pl = -1):
    """
    scattnlay(x, m[, theta, nmax, pl])

    Calculate the actual scattering parameters and amplitudes.

        x: Size parameters (2D ndarray)
        m: Relative refractive indices (2D ndarray)
        theta: Scattering angles where the scattering amplitudes will be
               calculated (optional, 1D ndarray)
        nmax: Maximum number of multipolar expansion terms to be used for the
              calculations. Only use it if you know what you are doing.
        pl: Index of PEC layer.

    Returns: (terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2)
    with
        terms: Number of multipolar expansion terms used for the calculations
        Qext: Efficiency factor for extinction
        Qsca: Efficiency factor for scattering
        Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)
        Qbk: Efficiency factor for backscattering
        Qpr: Efficiency factor for the radiation pressure
        g: Asymmetry factor (g = (Qext-Qpr)/Qsca)
        Albedo: Single scattering albedo (Albedo = Qsca/Qext)
        S1, S2: Complex scattering amplitudes
    """
    cdef Py_ssize_t i

    cdef np.ndarray[np.int_t, ndim = 1] terms = np.zeros(x.shape[0], dtype = np.int)

    cdef np.ndarray[np.float64_t, ndim = 1] Qext = np.zeros(x.shape[0], dtype = np.float64)
    cdef np.ndarray[np.float64_t, ndim = 1] Qabs = np.zeros(x.shape[0], dtype = np.float64)
    cdef np.ndarray[np.float64_t, ndim = 1] Qsca = np.zeros(x.shape[0], dtype = np.float64)
    cdef np.ndarray[np.float64_t, ndim = 1] Qbk = np.zeros(x.shape[0], dtype = np.float64)
    cdef np.ndarray[np.float64_t, ndim = 1] Qpr = np.zeros(x.shape[0], dtype = np.float64)
    cdef np.ndarray[np.float64_t, ndim = 1] g = np.zeros(x.shape[0], dtype = np.float64)
    cdef np.ndarray[np.float64_t, ndim = 1] Albedo = np.zeros(x.shape[0], dtype = np.float64)

    cdef np.ndarray[np.complex128_t, ndim = 2] S1 = np.zeros((x.shape[0], theta.shape[0]), dtype = np.complex128)
    cdef np.ndarray[np.complex128_t, ndim = 2] S2 = np.zeros((x.shape[0], theta.shape[0]), dtype = np.complex128)

    cdef np.ndarray[np.float64_t, ndim = 1] S1r
    cdef np.ndarray[np.float64_t, ndim = 1] S1i
    cdef np.ndarray[np.float64_t, ndim = 1] S2r
    cdef np.ndarray[np.float64_t, ndim = 1] S2i

    for i in range(x.shape[0]):
        S1r = np.zeros(theta.shape[0], dtype = np.float64)
        S1i = np.zeros(theta.shape[0], dtype = np.float64)
        S2r = np.zeros(theta.shape[0], dtype = np.float64)
        S2i = np.zeros(theta.shape[0], dtype = np.float64)

        terms[i] = nMie(x.shape[1], pl, x[i].copy('C'), m[i].copy('C'), theta.shape[0], theta.copy('C'), nmax, &Qext[i], &Qsca[i], &Qabs[i], &Qbk[i], &Qpr[i], &g[i], &Albedo[i], npy2c(S1r), npy2c(S1i), npy2c(S2r), npy2c(S2i))

        S1[i] = S1r.copy('C') + 1.0j*S1i.copy('C')
        S2[i] = S2r.copy('C') + 1.0j*S2i.copy('C')

    return terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2

def fieldnlay(np.ndarray[np.float64_t, ndim = 2] x, np.ndarray[np.complex128_t, ndim = 2] m, np.ndarray[np.float64_t, ndim = 2] coords, np.int_t nmax = -1, np.int_t pl = -1):
    """
    fieldnlay(x, m, coords[, theta, nmax, pl])

    Calculate the actual scattering parameters and amplitudes.

        x: Size parameters (2D ndarray)
        m: Relative refractive indices (2D ndarray)
        coords: Array containing all coordinates where the complex electric
                and magnetic fields will be calculated (2D ndarray)
        nmax: Maximum number of multipolar expansion terms to be used for the
              calculations. Only use it if you know what you are doing.
        pl: Index of PEC layer.

    Returns: (terms, E, H)
    with
        terms: Number of multipolar expansion terms used for the calculations
        E, H: Complex electric and magnetic field at the provided coordinates
    """
    cdef Py_ssize_t i

    cdef np.ndarray[np.int_t, ndim = 1] terms = np.zeros(x.shape[0], dtype = np.int)

    cdef np.ndarray[np.complex128_t, ndim = 3] E = np.zeros((x.shape[0], coords.shape[0], 3), dtype = np.complex128)
    cdef np.ndarray[np.complex128_t, ndim = 3] H = np.zeros((x.shape[0], coords.shape[0], 3), dtype = np.complex128)

    cdef np.ndarray[np.float64_t, ndim = 1] Erx
    cdef np.ndarray[np.float64_t, ndim = 1] Ery
    cdef np.ndarray[np.float64_t, ndim = 1] Erz
    cdef np.ndarray[np.float64_t, ndim = 1] Eix
    cdef np.ndarray[np.float64_t, ndim = 1] Eiy
    cdef np.ndarray[np.float64_t, ndim = 1] Eiz
    cdef np.ndarray[np.float64_t, ndim = 1] Hrx
    cdef np.ndarray[np.float64_t, ndim = 1] Hry
    cdef np.ndarray[np.float64_t, ndim = 1] Hrz
    cdef np.ndarray[np.float64_t, ndim = 1] Hix
    cdef np.ndarray[np.float64_t, ndim = 1] Hiy
    cdef np.ndarray[np.float64_t, ndim = 1] Hiz

    for i in range(x.shape[0]):
        Erx = np.zeros(coords.shape[0], dtype = np.float64)
        Ery = np.zeros(coords.shape[0], dtype = np.float64)
        Erz = np.zeros(coords.shape[0], dtype = np.float64)
        Eix = np.zeros(coords.shape[0], dtype = np.float64)
        Eiy = np.zeros(coords.shape[0], dtype = np.float64)
        Eiz = np.zeros(coords.shape[0], dtype = np.float64)
        Hrx = np.zeros(coords.shape[0], dtype = np.float64)
        Hry = np.zeros(coords.shape[0], dtype = np.float64)
        Hrz = np.zeros(coords.shape[0], dtype = np.float64)
        Hix = np.zeros(coords.shape[0], dtype = np.float64)
        Hiy = np.zeros(coords.shape[0], dtype = np.float64)
        Hiz = np.zeros(coords.shape[0], dtype = np.float64)

        terms[i] = nField(x.shape[1], pl, x[i].copy('C'), m[i].copy('C'), nmax, coords.shape[0], coords[:, 0].copy('C'), coords[:, 1].copy('C'), coords[:, 2].copy('C'), npy2c(Erx), npy2c(Ery), npy2c(Erz), npy2c(Eix), npy2c(Eiy), npy2c(Eiz), npy2c(Hrx), npy2c(Hry), npy2c(Hrz), npy2c(Hix), npy2c(Hiy), npy2c(Hiz))

        E[i] = np.vstack((Erx.copy('C') + 1.0j*Eix.copy('C'), Ery.copy('C') + 1.0j*Eiy.copy('C'), Erz.copy('C') + 1.0j*Eiz.copy('C'))).transpose()
        H[i] = np.vstack((Hrx.copy('C') + 1.0j*Hix.copy('C'), Hry.copy('C') + 1.0j*Hiy.copy('C'), Hrz.copy('C') + 1.0j*Hiz.copy('C'))).transpose()

    return terms, E, H

