#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2019 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#    Copyright (C) 2013-2019 Konstantin Ladutenko <kostyfisik@gmail.com>
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

import numpy as np

def scattcoeffs(x, m, nmax=-1, pl=-1, mp=False):
    """
    scattcoeffs(x, m[, nmax, pl, mp])

    Calculate the scattering coefficients required to calculate both the
    near- and far-field parameters.

        x: Size parameters (1D or 2D ndarray)
        m: Relative refractive indices (1D or 2D ndarray)
        nmax: Maximum number of multipolar expansion terms to be used for the
              calculations. Only use it if you know what you are doing, otherwise
              set this parameter to -1 and the function will calculate it.
        pl: Index of PEC layer. If there is none just send -1.
        mp: Use multiple (True) or double (False) precision.

    Returns: (terms, an, bn)
    with
        terms: Number of multipolar expansion terms used for the calculations
        an, bn: Complex scattering coefficients
    """

    if mp:
        from scattnlay_mp import scattcoeffs as scattcoeffs_
    else:
        from scattnlay_dp import scattcoeffs as scattcoeffs_

    if len(m.shape) != 1 and len(m.shape) != 2:
        raise ValueError('The relative refractive index (m) should be a 1-D or 2-D NumPy array.')
    if len(x.shape) == 1:
        if len(m.shape) == 1:
            return scattcoeffs_(x, m, nmax=nmax, pl=pl)
        else:
            raise ValueError('The number of of dimensions for the relative refractive index (m) and for the size parameter (x) must be equal.')
    elif len(x.shape) != 2:
        raise ValueError('The size parameter (x) should be a 1-D or 2-D NumPy array.')

    # Repeat the same m for all wavelengths
    if len(m.shape) == 1:
        m = np.repeat(m[np.newaxis, :], x.shape[0], axis=0)

    if nmax == -1:
        nstore = 0
    else:
        nstore = nmax

    terms = np.zeros((x.shape[0]), dtype=int)
    an = np.zeros((0, nstore), dtype=complex)
    bn = np.zeros((0, nstore), dtype=complex)

    for i, xi in enumerate(x):
        terms[i], a, b = scattcoeffs_(xi, m[i], nmax=nmax, pl=pl)

        if terms[i] > nstore:
            nstore = terms[i]
            an.resize((an.shape[0], nstore))
            bn.resize((bn.shape[0], nstore))

        an = np.vstack((an, a))
        bn = np.vstack((bn, b))

    return terms, an, bn
#scattcoeffs()

def expancoeffs(x, m, nmax=-1, pl=-1, mp=False):
    """
    expancoeffs(x, m[, nmax, pl, mp])

    Calculate the scattering coefficients required to calculate both the
    near- and far-field parameters.

        x: Size parameters (1D or 2D ndarray)
        m: Relative refractive indices (1D or 2D ndarray)
        nmax: Maximum number of multipolar expansion terms to be used for the
              calculations. Only use it if you know what you are doing, otherwise
              set this parameter to -1 and the function will calculate it.
        pl: Index of PEC layer. If there is none just send -1.
        mp: Use multiple (True) or double (False) precision.

    Returns: (terms, an, bn, cn, dn)
    with
        terms: Number of multipolar expansion terms used for the calculations
        an, bn, cn, dn: Complex expansion coefficients of each layer
    """

    if mp:
        from scattnlay_mp import expancoeffs as expancoeffs_
    else:
        from scattnlay_dp import expancoeffs as expancoeffs_

    if len(m.shape) != 1 and len(m.shape) != 2:
        raise ValueError('The relative refractive index (m) should be a 1-D or 2-D NumPy array.')
    if len(x.shape) == 1:
        if len(m.shape) == 1:
            return expancoeffs_(x, m, nmax=nmax, pl=pl)
        else:
            raise ValueError('The number of of dimensions for the relative refractive index (m) and for the size parameter (x) must be equal.')
    elif len(x.shape) != 2:
        raise ValueError('The size parameter (x) should be a 1-D or 2-D NumPy array.')

    # Repeat the same m for all wavelengths
    if len(m.shape) == 1:
        m = np.repeat(m[np.newaxis, :], x.shape[0], axis=0)

    if nmax == -1:
        nstore = 0
    else:
        nstore = nmax

    terms = np.zeros((x.shape[0]), dtype=int)
    an = np.zeros((0, x.shape[1], nstore), dtype=complex)
    bn = np.zeros((0, x.shape[1], nstore), dtype=complex)
    cn = np.zeros((0, x.shape[1], nstore), dtype=complex)
    dn = np.zeros((0, x.shape[1], nstore), dtype=complex)

    for i, xi in enumerate(x):
        terms[i], a, b, c, d = expancoeffs_(xi, m[i], nmax=nmax, pl=pl)

        if terms[i] > nstore:
            nstore = terms[i]
            an.resize((an.shape[0], an.shape[1], nstore))
            bn.resize((bn.shape[0], bn.shape[1], nstore))
            cn.resize((cn.shape[0], cn.shape[1], nstore))
            dn.resize((dn.shape[0], dn.shape[1], nstore))

        an = np.vstack((an, a))
        bn = np.vstack((bn, b))
        cn = np.vstack((cn, c))
        dn = np.vstack((dn, d))

    return terms, an, bn, cn, dn
#expancoeffs()


def scattnlay(x, m, theta=np.zeros(0, dtype=float), nmax=-1, pl=-1, mp=False):
    """
    scattnlay(x, m[, theta, nmax, pl, mp])

    Calculate the actual scattering parameters and amplitudes.

        x: Size parameters (1D or 2D ndarray)
        m: Relative refractive indices (1D or 2D ndarray)
        theta: Scattering angles where the scattering amplitudes will be
               calculated (optional, 1D ndarray)
        nmax: Maximum number of multipolar expansion terms to be used for the
              calculations. Only use it if you know what you are doing.
        pl: Index of PEC layer. If there is none just send -1.
        mp: Use multiple (True) or double (False) precision.

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

    if mp:
        from scattnlay_mp import scattnlay as scattnlay_
    else:
        from scattnlay_dp import scattnlay as scattnlay_

    if len(m.shape) != 1 and len(m.shape) != 2:
        raise ValueError('The relative refractive index (m) should be a 1-D or 2-D NumPy array.')
    if len(x.shape) == 1:
        if len(m.shape) == 1:
            return scattnlay_(x, m, theta, nmax=nmax, pl=pl)
        else:
            raise ValueError('The number of of dimensions for the relative refractive index (m) and for the size parameter (x) must be equal.')
    elif len(x.shape) != 2:
        raise ValueError('The size parameter (x) should be a 1-D or 2-D NumPy array.')
    if len(theta.shape) != 1:
        raise ValueError('The scattering angles (theta) should be a 1-D NumPy array.')

    # Repeat the same m for all wavelengths
    if len(m.shape) == 1:
        m = np.repeat(m[np.newaxis, :], x.shape[0], axis=0)

    terms = np.zeros((x.shape[0]), dtype=int)
    Qext = np.zeros((x.shape[0]), dtype=float)
    Qsca = np.zeros((x.shape[0]), dtype=float)
    Qabs = np.zeros((x.shape[0]), dtype=float)
    Qbk = np.zeros((x.shape[0]), dtype=float)
    Qpr = np.zeros((x.shape[0]), dtype=float)
    g = np.zeros((x.shape[0]), dtype=float)
    Albedo = np.zeros((x.shape[0]), dtype=float)
    S1 = np.zeros((x.shape[0], theta.shape[0]), dtype=complex)
    S2 = np.zeros((x.shape[0], theta.shape[0]), dtype=complex)

    for i, xi in enumerate(x):
        terms[i], Qext[i], Qsca[i], Qabs[i], Qbk[i], Qpr[i], g[i], Albedo[i], S1[i], S2[i] = scattnlay_(xi, m[i], theta, nmax=nmax, pl=pl)

    return terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2
#scattnlay()


def fieldnlay(x, m, xp, yp, zp, nmax=-1, pl=-1, mp=False):
    """
    fieldnlay(x, m, xp, yp, zp[, nmax, pl, mp])

    Calculate the actual scattering parameters and amplitudes.

        x: Size parameters (1D or 2D ndarray)
        m: Relative refractive indices (1D or 2D ndarray)
        xp: Array containing all X coordinates to calculate the complex
            electric and magnetic fields (1D* ndarray)
        yp: Array containing all Y coordinates to calculate the complex
            electric and magnetic fields (1D* ndarray)
        zp: Array containing all Z coordinates to calculate the complex
            electric and magnetic fields (1D* ndarray)
        nmax: Maximum number of multipolar expansion terms to be used for the
              calculations. Only use it if you know what you are doing.
        pl: Index of PEC layer. If there is none just send -1.
        mp: Use multiple (True) or double (False) precision.

    Returns: (terms, E, H)
    with
        terms: Number of multipolar expansion terms used for the calculations
        E, H: Complex electric and magnetic field at the provided coordinates

    *Note: We assume that the coordinates are referred to the first wavelength
           (or structure) and correct it for the following ones
    """

    if mp:
        from scattnlay_mp import fieldnlay as fieldnlay_
    else:
        from scattnlay_dp import fieldnlay as fieldnlay_

    if len(m.shape) != 1 and len(m.shape) != 2:
        raise ValueError('The relative refractive index (m) should be a 1-D or 2-D NumPy array.')
    if len(x.shape) == 1:
        if len(m.shape) == 1:
            return fieldnlay_(x, m, xp, yp, zp, nmax=nmax, pl=pl)
        else:
            raise ValueError('The number of of dimensions for the relative refractive index (m) and for the size parameter (x) must be equal.')
    elif len(x.shape) != 2:
        raise ValueError('The size parameter (x) should be a 1-D or 2-D NumPy array.')

    # Repeat the same m for all wavelengths
    if len(m.shape) == 1:
        m = np.repeat(m[np.newaxis, :], x.shape[0], axis=0)

    terms = np.zeros((x.shape[0]), dtype=int)
    E = np.zeros((x.shape[0], xp.shape[0], 3), dtype=complex)
    H = np.zeros((x.shape[0], xp.shape[0], 3), dtype=complex)

    for i, xi in enumerate(x):
        # (2020/05/12) We assume that the coordinates are referred to the first wavelength
        #              (or structure) and correct it for the following ones
        terms[i], E[i], H[i] = fieldnlay_(xi, m[i], xp*xi[-1]/x[0, -1], yp*xi[-1]/x[0, -1], zp*xi[-1]/x[0, -1], nmax=nmax, pl=pl)

    return terms, E, H
#fieldnlay()

