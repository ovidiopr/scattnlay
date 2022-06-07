#ifndef SRC_NMIE_BASIC_HPP_
#define SRC_NMIE_BASIC_HPP_
//******************************************************************************
//    Copyright (C) 2009-2018  Ovidio Pena <ovidio@bytesfall.com>
//    Copyright (C) 2013-2018  Konstantin Ladutenko <kostyfisik@gmail.com>
//
//    This file is part of scattnlay
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    The only additional remark is that we expect that all publications
//    describing work using this software, or all commercial products
//    using it, cite at least one of the following references:
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by
//        a multilayered sphere," Computer Physics Communications,
//        vol. 180, Nov. 2009, pp. 2348-2354.
//    [2] K. Ladutenko, U. Pal, A. Rivera, and O. Pena-Rodriguez, "Mie
//        calculation of electromagnetic near-field for a multilayered
//        sphere," Computer Physics Communications, vol. 214, May 2017,
//        pp. 225-230.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//******************************************************************************

//******************************************************************************
// This class implements the algorithm for a multilayered sphere described by:
//    [1] W. Yang, "Improved recursive algorithm for light scattering by a
//        multilayered sphere,” Applied Optics, vol. 42, Mar. 2003, pp. 1710-1720.
//
// You can find the description of all the used equations in:
//    [2] O. Pena and U. Pal, "Scattering of electromagnetic radiation by
//        a multilayered sphere," Computer Physics Communications,
//        vol. 180, Nov. 2009, pp. 2348-2354.
//    [3] K. Ladutenko, U. Pal, A. Rivera, and O. Pena-Rodriguez, "Mie
//        calculation of electromagnetic near-field for a multilayered
//        sphere," Computer Physics Communications, vol. 214, May 2017,
//        pp. 225-230.
//
// Hereinafter all equations numbers refer to [2]
//******************************************************************************
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "nmie.hpp"
#include "special-functions-impl.hpp"

namespace nmie {
template <typename FloatType>
void MesoMie<FloatType>::calc_ab(int nmax,
                                 FloatType R,
                                 FloatType xd,
                                 FloatType xm,
                                 FloatType eps_d,
                                 FloatType eps_m,
                                 FloatType d_parallel,
                                 FloatType d_perp) {
  an_.resize(nmax, static_cast<FloatType>(0.0));
  bn_.resize(nmax, static_cast<FloatType>(0.0));

  std::vector<std::complex<FloatType>> D1_xd(nmax + 1), D1_xm(nmax + 1),
      D3_xd(nmax + 1), D3_xm(nmax + 1);
  for (int n = 0; n <= nmax_; n++) {
    D1_xd[n] = std::complex<FloatType>(0.0, -1.0);
    D1_xm[n] = std::complex<FloatType>(0.0, -1.0);
    D3_xd[n] = std::complex<FloatType>(0.0, 1.0);
  }

  std::vector<std::complex<FloatType>> Psi_xd(nmax + 1), Zeta_xd(nmax + 1),
      Psi_xm(nmax + 1), Zeta_xm(nmax + 1), PsiZeta_xd(nmax + 1),
      PsiZeta_xm(nmax + 1);

  std::complex<FloatType> cxd(xd, 0), cxm(xm, 0);
  evalDownwardD1(cxd, D1_xd);
  evalUpwardPsi(cxd, D1_xd, Psi_xd);
  evalUpwardD3(cxd, D1_xd, D3_xd, PsiZeta_xd);
  for (unsigned int i = 0; i < Zeta.size(); i++) {
    Zeta[i] = PsiZeta[i] / Psi[i];
  }

  evalDownwardD1(cxm, D1_xm);
  evalUpwardD3(cxm, D1_xm, D3_xm, PsiZeta_xm);
}

}  // namespace nmie
#endif
