#    Copyright (C) 2009-2015  Ovidio Pena <ovidio@bytesfall.com>
#
#    This file is part of python-scattnlay
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
#    using it, cite the following reference:
#    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by
#        a multilayered sphere," Computer Physics Communications,
#        vol. 180, Nov. 2009, pp. 2348-2354.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division
import numpy as np
cimport numpy as np

cdef extern from "<vector>" namespace "std":
    cdef cppclass vector[T]:
        cppclass iterator:
            T operator*()
            iterator operator++()
            bint operator==(iterator)
            bint operator!=(iterator)
        vector()
        void push_back(T&)
        T& operator[](int)
        T& at(int)
        iterator begin()
        iterator end()

cdef inline double *npy2c(np.ndarray a):
    assert a.dtype == np.float64

    if not (<object>a).flags["C_CONTIGUOUS"]: # Array is not contiguous, need to make contiguous copy
        a = a.copy('C')

    # Return data pointer
    return <double *>(a.data)

cdef extern from "py_nmie.h":
    cdef int nMie(int L, int pl, vector[double] x, vector[complex] m, int nTheta, vector[double] Theta, int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, double S1r[], double S1i[], double S2r[], double S2i[])
    cdef int nField(int L, int pl, vector[double] x, vector[complex] m, int nmax, int nCoords, vector[double] Xp, vector[double] Yp, vector[double] Zp, double Er[], double Ei[], double Hr[], double Hi[])

def scattnlay(np.ndarray[np.float64_t, ndim = 2] x, np.ndarray[np.complex128_t, ndim = 2] m, np.ndarray[np.float64_t, ndim = 1] theta = np.zeros(0, dtype = np.float64), np.int_t pl = -1, np.int_t nmax = -1):
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

def fieldnlay(np.ndarray[np.float64_t, ndim = 2] x, np.ndarray[np.complex128_t, ndim = 2] m, np.ndarray[np.float64_t, ndim = 2] coords = np.zeros((0, 3), dtype = np.float64), np.int_t pl = 0, np.int_t nmax = 0):
    cdef Py_ssize_t i

    cdef np.ndarray[np.int_t, ndim = 1] terms = np.zeros(x.shape[0], dtype = np.int)

    cdef np.ndarray[np.complex128_t, ndim = 2] E = np.zeros((x.shape[0], coords.shape[0]), dtype = np.complex128)
    cdef np.ndarray[np.complex128_t, ndim = 2] H = np.zeros((x.shape[0], coords.shape[0]), dtype = np.complex128)

    cdef np.ndarray[np.float64_t, ndim = 1] Er
    cdef np.ndarray[np.float64_t, ndim = 1] Ei
    cdef np.ndarray[np.float64_t, ndim = 1] Hr
    cdef np.ndarray[np.float64_t, ndim = 1] Hi

    for i in range(x.shape[0]):
        Er = np.zeros(coords.shape[0], dtype = np.float64)
        Ei = np.zeros(coords.shape[0], dtype = np.float64)
        Hr = np.zeros(coords.shape[0], dtype = np.float64)
        Hi = np.zeros(coords.shape[0], dtype = np.float64)

        terms[i] = nField(x.shape[1], pl, x[i].copy('C'), m[i].copy('C'), nmax, coords.shape[0], coords[:, 0].copy('C'), coords[:, 1].copy('C'), coords[:, 2].copy('C'), npy2c(Er), npy2c(Ei), npy2c(Hr), npy2c(Hi))

        E[i] = Er.copy('C') + 1.0j*Ei.copy('C')
        H[i] = Hr.copy('C') + 1.0j*Hi.copy('C')

    return terms, E, H

