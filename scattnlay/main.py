#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2021 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#    Copyright (C) 2013-2021 Konstantin Ladutenko <kostyfisik@gmail.com>
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
try:
    from .scattnlay_dp import mie_dp, mesomie_dp
except ImportError:
    from scattnlay_dp import mie_dp, mesomie_dp
import numpy as np
import sys

mie_mp = None
# mesomie_mp = None
try:
    try:
        from .scattnlay_mp import mie_mp as mie_mp_
    except ImportError:
        from scattnlay_mp import mie_mp as mie_mp_
    mie_mp = mie_mp_()
    # from scattnlay_mp import mesomie_mp as mesomie_mp_
    # mesomie_mp = mesomie_mp_()
except:
    pass
mie_simd = None
try:
    try:
        from . import scattnlay_simd
    except ImportError:
        import scattnlay_simd
    # Create a compatibility wrapper so it behaves like mie_dp/mie_mp classes
    class MieSIMDWrapper:
        def RunMieBatch(self, x, m, theta=None):
            if theta is None:
                return scattnlay_simd.RunMieBatch(np.atleast_1d(x), np.atleast_1d(m))
            else:
                return scattnlay_simd.RunMieBatch(np.atleast_1d(x), np.atleast_1d(m), np.atleast_1d(theta))
        
        # Add basic methods for test parity
        def SetLayersSize(self, x): self._x = x
        def SetLayersIndex(self, m): self._m = m
        def SetAngles(self, theta): self._theta = theta
        def RunMieCalculation(self):
            # For bulk test parity, simulate class behavior
            if hasattr(self, '_theta'):
                self._res = scattnlay_simd.RunMieBatch(np.atleast_1d(self._x), np.atleast_1d(self._m), np.atleast_1d(self._theta))
            else:
                self._res = scattnlay_simd.RunMieBatch(np.atleast_1d(self._x), np.atleast_1d(self._m))
        def GetQext(self): return self._res['Qext'][0]
        def GetQsca(self): return self._res['Qsca'][0]
        def GetS1(self): return self._res['S1'][0]
        def GetS2(self): return self._res['S2'][0]

    mie_simd = MieSIMDWrapper()
except ImportError:
    pass

mie = mie_dp()
mesomie = mesomie_dp()


def scattcoeffs_(x, m, nmax=-1, pl=-1, mp=False):
    if mp and mie_mp:
        from scattnlay_mp import mie_mp as mie_
    else:
        if mp:
            print('Failed to load multiprecision module, using double precision instead...',
                  file=sys.stderr)
        from scattnlay_dp import mie_dp as mie_
        # from scattnlay_mp import mie_mp as mie_
    mie = mie_()
    mie.SetLayersSize(x)
    mie.SetLayersIndex(m)
    mie.SetPECLayer(pl)
    mie.SetMaxTerms(nmax)
    mie.calcScattCoeffs()
    terms = mie.GetMaxTerms()
    a = mie.GetAn()
    b = mie.GetBn()
    return terms, a, b


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

    if len(m.shape) != 1 and len(m.shape) != 2:
        raise ValueError(
            'The relative refractive index (m) should be a 1-D or 2-D NumPy array.')
    if len(x.shape) == 1:
        if len(m.shape) == 1:
            return scattcoeffs_(x, m, nmax=nmax, pl=pl, mp=mp)
        else:
            raise ValueError(
                'The number of of dimensions for the relative refractive index (m) and for the size parameter (x) must be equal.')
    elif len(x.shape) != 2:
        raise ValueError(
            'The size parameter (x) should be a 1-D or 2-D NumPy array.')

    # Repeat the same m for all wavelengths
    if len(m.shape) == 1:
        m = np.repeat(m[np.newaxis, :], x.shape[0], axis=0)

    if nmax == -1:
        if x.size > 0:
            max_x = np.max(x)
            nstore = int(max_x + 4.0 * max_x**(1.0/3.0) + 2.0) + 10
        else:
            nstore = 0
    else:
        nstore = nmax

    num_points = x.shape[0]
    terms = np.zeros(num_points, dtype=int)
    an = np.zeros((num_points, nstore), dtype=complex)
    bn = np.zeros((num_points, nstore), dtype=complex)

    for i, xi in enumerate(x):
        terms[i], a, b = scattcoeffs_(xi, m[i], nmax=nmax, pl=pl, mp=mp)

        if terms[i] > nstore:
            nstore = terms[i]
            an.resize((num_points, nstore), refcheck=False)
            bn.resize((num_points, nstore), refcheck=False)

        an[i, :terms[i]] = a
        bn[i, :terms[i]] = b

    return terms, an, bn


def expancoeffs_(x, m, nmax=-1, pl=-1, mp=False):
    if mp and mie_mp:
        from scattnlay_mp import mie_mp as mie_
    else:
        if mp:
            print('Failed to load multiprecision module, using double precision instead...',
                  file=sys.stderr)
        from scattnlay_dp import mie_dp as mie_
        # from scattnlay_mp import mie_mp as mie_
    mie = mie_()
    mie.SetLayersSize(x)
    mie.SetLayersIndex(m)
    mie.SetPECLayer(pl)
    mie.SetMaxTerms(nmax)
    mie.calcScattCoeffs()
    mie.calcExpanCoeffs()
    terms = mie.GetMaxTerms()
    an = mie.GetLayerAn()
    bn = mie.GetLayerBn()
    cn = mie.GetLayerCn()
    dn = mie.GetLayerDn()
    return terms, an, bn, cn, dn


# TODO verify that expancoeffs() is really working
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
    if len(m.shape) != 1 and len(m.shape) != 2:
        raise ValueError(
            'The relative refractive index (m) should be a 1-D or 2-D NumPy array.')
    if len(x.shape) == 1:
        if len(m.shape) == 1:
            return expancoeffs_(x, m, nmax=nmax, pl=pl, mp=mp)
        else:
            raise ValueError(
                'The number of of dimensions for the relative refractive index (m) and for the size parameter (x) must be equal.')
    elif len(x.shape) != 2:
        raise ValueError(
            'The size parameter (x) should be a 1-D or 2-D NumPy array.')

    # Repeat the same m for all wavelengths
    if len(m.shape) == 1:
        m = np.repeat(m[np.newaxis, :], x.shape[0], axis=0)

    if nmax == -1:
        if x.size > 0:
            max_x = np.max(x)
            nstore = int(max_x + 4.0 * max_x**(1.0/3.0) + 2.0) + 10
        else:
            nstore = 0
    else:
        nstore = nmax

    num_points = x.shape[0]
    num_layers = x.shape[1]
    terms = np.zeros(num_points, dtype=int)
    an = np.zeros((num_points, num_layers + 1, nstore), dtype=complex)
    bn = np.zeros((num_points, num_layers + 1, nstore), dtype=complex)
    cn = np.zeros((num_points, num_layers + 1, nstore), dtype=complex)
    dn = np.zeros((num_points, num_layers + 1, nstore), dtype=complex)

    for i, xi in enumerate(x):
        terms[i], a, b, c, d = expancoeffs_(xi, m[i], nmax=nmax, pl=pl, mp=mp)

        if terms[i] > nstore:
            nstore = terms[i]
            an.resize((num_points, num_layers + 1, nstore), refcheck=False)
            bn.resize((num_points, num_layers + 1, nstore), refcheck=False)
            cn.resize((num_points, num_layers + 1, nstore), refcheck=False)
            dn.resize((num_points, num_layers + 1, nstore), refcheck=False)

        an[i, :, :terms[i]] = a
        bn[i, :, :terms[i]] = b
        cn[i, :, :terms[i]] = c
        dn[i, :, :terms[i]] = d

    return terms, an, bn, cn, dn


def scattnlay_(x, m, theta=np.zeros(0, dtype=float), nmax=-1, pl=-1, mp=False):
    if mp and mie_mp:
        from scattnlay_mp import mie_mp as mie_
    else:
        if mp:
            print('Failed to load multiprecision module, using double precision instead...',
                  file=sys.stderr)
        from scattnlay_dp import mie_dp as mie_
    mie = mie_()
    mie.SetLayersSize(x)
    mie.SetLayersIndex(m)
    mie.SetAngles(theta)
    mie.SetPECLayer(pl)
    mie.SetMaxTerms(nmax)
    mie.RunMieCalculation()
    Qext = mie.GetQext()
    Qsca = mie.GetQsca()
    Qabs = mie.GetQabs()
    Qbk = mie.GetQbk()
    Qpr = mie.GetQpr()
    g = mie.GetAsymmetryFactor()
    Albedo = mie.GetAlbedo()
    terms = mie.GetMaxTerms()
    S1 = mie.GetS1()
    S2 = mie.GetS2()
    return terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2


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
    if len(m.shape) != 1 and len(m.shape) != 2:
        raise ValueError(
            'The relative refractive index (m) should be a 1-D or 2-D NumPy array.')
    if len(x.shape) == 1:
        if len(m.shape) == 1:
            return scattnlay_(x, m, theta, nmax=nmax, pl=pl, mp=mp)
        else:
            raise ValueError(
                'The number of of dimensions for the relative refractive index (m) and for the size parameter (x) must be equal.')
    elif len(x.shape) != 2:
        raise ValueError(
            'The size parameter (x) should be a 1-D or 2-D NumPy array.')
    if len(theta.shape) != 1:
        raise ValueError(
            'The scattering angles (theta) should be a 1-D NumPy array.')

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
        terms[i], Qext[i], Qsca[i], Qabs[i], Qbk[i], Qpr[i], g[i], Albedo[i], S1[i], S2[i] = scattnlay_(
            xi, m[i], theta, nmax=nmax, pl=pl, mp=mp)

    return terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2


def fieldnlay_(x, m, xp, yp, zp, nmax=-1, pl=-1, mp=False):
    if mp and mie_mp:
        from scattnlay_mp import mie_mp as mie_
    else:
        if mp:
            print('Failed to load multiprecision module, using double precision instead...',
                  file=sys.stderr)
        from scattnlay_dp import mie_dp as mie_
        # from scattnlay_mp import mie_mp as mie_
    mie = mie_()
    mie.SetLayersSize(x)
    mie.SetLayersIndex(m)
    mie.SetPECLayer(pl)
    mie.SetMaxTerms(nmax)
    mie.SetFieldCoords(xp, yp, zp)
    mie.RunFieldCalculation()
    terms = mie.GetMaxTerms()
    E = mie.GetFieldE()
    H = mie.GetFieldH()
    return terms, E, H


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
    if len(m.shape) != 1 and len(m.shape) != 2:
        raise ValueError(
            'The relative refractive index (m) should be a 1-D or 2-D NumPy array.')
    if len(x.shape) == 1:
        if len(m.shape) == 1:
            return fieldnlay_(x, m, xp, yp, zp, nmax=nmax, pl=pl, mp=mp)
        else:
            raise ValueError(
                'The number of of dimensions for the relative refractive index (m) and for the size parameter (x) must be equal.')
    elif len(x.shape) != 2:
        raise ValueError(
            'The size parameter (x) should be a 1-D or 2-D NumPy array.')

    # Repeat the same m for all wavelengths
    if len(m.shape) == 1:
        m = np.repeat(m[np.newaxis, :], x.shape[0], axis=0)

    terms = np.zeros((x.shape[0]), dtype=int)
    E = np.zeros((x.shape[0], xp.shape[0], 3), dtype=complex)
    H = np.zeros((x.shape[0], xp.shape[0], 3), dtype=complex)

    for i, xi in enumerate(x):
        # (2020/05/12) We assume that the coordinates are referred to the first wavelength
        #              (or structure) and correct it for the following ones
        terms[i], E[i], H[i] = fieldnlay_(
            xi, m[i], xp*xi[-1]/x[0, -1], yp*xi[-1]/x[0, -1], zp*xi[-1]/x[0, -1], nmax=nmax, pl=pl, mp=mp)

    return terms, E, H
