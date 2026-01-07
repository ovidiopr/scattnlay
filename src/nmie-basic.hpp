#ifndef SRC_NMIE_BASIC_HPP_
#define SRC_NMIE_BASIC_HPP_
//***************************************************************************//
//    Copyright (C) 2009-2022  Ovidio Pena <ovidio@bytesfall.com>            //
//    Copyright (C) 2013-202  Konstantin Ladutenko <kostyfisik@gmail.com>    //
//                                                                           //
//    This file is part of scattnlay                                         //
//                                                                           //
//    This program is free software: you can redistribute it and/or modify   //
//    it under the terms of the GNU General Public License as published by   //
//    the Free Software Foundation, either version 3 of the License, or      //
//    (at your option) any later version.                                    //
//                                                                           //
//    This program is distributed in the hope that it will be useful,        //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of         //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          //
//    GNU General Public License for more details.                           //
//                                                                           //
//    The only additional remark is that we expect that all publications     //
//    describing work using this software, or all commercial products        //
//    using it, cite at least one of the following references:               //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by    //
//        a multilayered sphere," Computer Physics Communications,           //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                //
//    [2] K. Ladutenko, U. Pal, A. Rivera, and O. Pena-Rodriguez, "Mie       //
//        calculation of electromagnetic near-field for a multilayered       //
//        sphere," Computer Physics Communications, vol. 214, May 2017,      //
//        pp. 225-230.                                                       //
//                                                                           //
//    You should have received a copy of the GNU General Public License      //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.  //
//***************************************************************************//
//***************************************************************************//
// This class implements the algorithm for a multilayered sphere described
//    by:
//    [1] W. Yang, "Improved recursive algorithm for light scattering by a
//        multilayered sphere,” Applied Optics, vol. 42, Mar. 2003, pp.
//        1710-1720.
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
//*****************************************************************************//
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "nmie.hpp"
#include "special-functions-impl.hpp"

namespace nmie {

template <typename FloatType, typename Engine = ScalarEngine<FloatType>, typename ComplexType>
ComplexType calc_an(typename Engine::RealV n_real,
                      typename Engine::RealV XL, 
                      ComplexType Ha,
                      ComplexType mL,
                      ComplexType PsiXL,
                      ComplexType ZetaXL,
                      ComplexType PsiXLM1,
                      ComplexType ZetaXLM1) {
  auto term1 = (Ha / mL) + (n_real / XL);
  auto Num = term1 * PsiXL - PsiXLM1;
  auto Denom = term1 * ZetaXL - ZetaXLM1;
  return Num / Denom;
}

template <typename FloatType, typename Engine = ScalarEngine<FloatType>, typename ComplexType>
ComplexType calc_an(int n,
                      typename Engine::RealV XL, 
                      ComplexType Ha,
                      ComplexType mL,
                      ComplexType PsiXL,
                      ComplexType ZetaXL,
                      ComplexType PsiXLM1,
                      ComplexType ZetaXLM1) {
    auto n_val = Engine::set(static_cast<FloatType>(n));
    return calc_an<FloatType, Engine, ComplexType>(n_val, XL, Ha, mL, PsiXL, ZetaXL, PsiXLM1, ZetaXLM1);
}

template <typename FloatType, typename Engine = ScalarEngine<FloatType>, typename ComplexType>
ComplexType calc_bn(typename Engine::RealV n_real,
                      typename Engine::RealV XL, 
                      ComplexType Hb,
                      ComplexType mL,
                      ComplexType PsiXL,
                      ComplexType ZetaXL,
                      ComplexType PsiXLM1,
                      ComplexType ZetaXLM1) {
  auto term1 = (mL * Hb) + (n_real / XL);
  auto Num = term1 * PsiXL - PsiXLM1;
  auto Denom = term1 * ZetaXL - ZetaXLM1;
  return Num / Denom;
}

template <typename FloatType, typename Engine = ScalarEngine<FloatType>, typename ComplexType>
ComplexType calc_bn(int n,
                      typename Engine::RealV XL, 
                      ComplexType Hb,
                      ComplexType mL,
                      ComplexType PsiXL,
                      ComplexType ZetaXL,
                      ComplexType PsiXLM1,
                      ComplexType ZetaXLM1) {
    auto n_val = Engine::set(static_cast<FloatType>(n));
    return calc_bn<FloatType, Engine, ComplexType>(n_val, XL, Hb, mL, PsiXL, ZetaXL, PsiXLM1, ZetaXLM1);
}

template <typename FloatType,
          typename Engine = ScalarEngine<FloatType>,
          typename ComplexType>
void computeAnBnBatch(typename Engine::RealV n_real,
                      typename Engine::RealV XL,
                      ComplexType Ha,
                      ComplexType Hb,
                      ComplexType mL,
                      ComplexType PsiXL,
                      ComplexType ZetaXL,
                      ComplexType PsiXLM1,
                      ComplexType ZetaXLM1,
                      ComplexType& an,
                      ComplexType& bn) {
  an = calc_an<FloatType, Engine, ComplexType>(n_real, XL, Ha, mL, PsiXL,
                                               ZetaXL, PsiXLM1, ZetaXLM1);
  bn = calc_bn<FloatType, Engine, ComplexType>(n_real, XL, Hb, mL, PsiXL,
                                               ZetaXL, PsiXLM1, ZetaXLM1);
}

template <typename FloatType, typename Engine = ScalarEngine<FloatType>>
void computeLayerCoeffsHelper(int nmax,
                              int l,
                              int pl,
                              typename Engine::RealV x_l,
                              typename Engine::RealV x_lm1,
                              typename Engine::ComplexV m_l,
                              typename Engine::ComplexV m_lm1,
                              typename Engine::ComplexV z1,
                              typename Engine::ComplexV z2,
                              const std::complex<FloatType>* D1_mlxl,
                              const std::complex<FloatType>* D3_mlxl,
                              const std::complex<FloatType>* D1_mlxlM1,
                              const std::complex<FloatType>* D3_mlxlM1,
                              std::complex<FloatType>* Q_l,
                              std::complex<FloatType>* Ha_l,
                              std::complex<FloatType>* Hb_l,
                              const std::complex<FloatType>* Ha_lm1,
                              const std::complex<FloatType>* Hb_lm1) {
  // Q[l][0] calculation
  auto z1_re = Engine::get_real(z1);
  auto z1_im = Engine::get_imag(z1);
  auto z2_re = Engine::get_real(z2);
  auto z2_im = Engine::get_imag(z2);

  auto minus_two = Engine::set(-2.0);
  auto exp_term_num = Engine::exp(minus_two * (z1_im - z2_im));

  auto arg_z2 = minus_two * z2_re;
  auto exp_z2 = Engine::exp(minus_two * z2_im);
  auto cos_z2 = Engine::cos(arg_z2);
  auto sin_z2 = Engine::sin(arg_z2);

  auto num_re = exp_term_num * (cos_z2 - exp_z2);
  auto num_im = exp_term_num * sin_z2;
  auto Num = Engine::make_complex(num_re, num_im);

  auto arg_z1 = minus_two * z1_re;
  auto exp_z1 = Engine::exp(minus_two * z1_im);
  auto cos_z1 = Engine::cos(arg_z1);
  auto sin_z1 = Engine::sin(arg_z1);

  auto denom_re = cos_z1 - exp_z1;
  auto denom_im = sin_z1;
  auto Denom = Engine::make_complex(denom_re, denom_im);

  auto Q_curr = Num / Denom;
  Engine::store(Q_curr, &Q_l[0]);

  auto ratio = x_lm1 / x_l;
  auto ratio_sq = ratio * ratio;
  auto ratio_sq_c = Engine::make_complex(ratio_sq, Engine::set(0.0));

  const size_t lanes = Engine::Lanes();

  for (int n = 1; n <= nmax; ++n) {
    auto n_val = Engine::set(static_cast<FloatType>(n));
    auto n_c = Engine::make_complex(n_val, Engine::set(0.0));

    auto d1_n = Engine::load(&D1_mlxl[n * lanes]);
    auto d3_nm1 = Engine::load(&D3_mlxl[(n - 1) * lanes]);

    auto term1 = (z1 * d1_n) + n_c;
    auto term2 = n_c - (z1 * d3_nm1);
    auto Num_n = term1 * term2;

    auto d1_m1_n = Engine::load(&D1_mlxlM1[n * lanes]);
    auto d3_m1_nm1 = Engine::load(&D3_mlxlM1[(n - 1) * lanes]);

    auto term3 = (z2 * d1_m1_n) + n_c;
    auto term4 = n_c - (z2 * d3_m1_nm1);
    auto Denom_n = term3 * term4;

    auto factor = ratio_sq_c * Q_curr;
    Q_curr = (factor * Num_n) / Denom_n;
    Engine::store(Q_curr, &Q_l[n * lanes]);

    // Ha
    typename Engine::ComplexV G1_ha, G2_ha;
    auto ha_prev = Engine::load(&Ha_lm1[(n - 1) * lanes]);

    if ((l - 1) == pl) {
      auto neg_one = Engine::set(-1.0);
      auto neg_one_c = Engine::make_complex(neg_one, Engine::set(0.0));
      G1_ha = d1_m1_n * neg_one_c;
      auto d3_m1_n = Engine::load(&D3_mlxlM1[n * lanes]);
      G2_ha = d3_m1_n * neg_one_c;
    } else {
      auto term_ha = m_l * ha_prev;
      G1_ha = term_ha - (m_lm1 * d1_m1_n);
      auto d3_m1_n = Engine::load(&D3_mlxlM1[n * lanes]);
      G2_ha = term_ha - (m_lm1 * d3_m1_n);
    }

    auto Temp_ha = Q_curr * G1_ha;
    auto d1_n_curr = Engine::load(&D1_mlxl[n * lanes]);
    auto d3_n = Engine::load(&D3_mlxl[n * lanes]);
    auto Num_ha = (G2_ha * d1_n_curr) - (Temp_ha * d3_n);
    auto Denom_ha = G2_ha - Temp_ha;
    auto Ha_curr = Num_ha / Denom_ha;
    Engine::store(Ha_curr, &Ha_l[(n - 1) * lanes]);

    // Hb
    typename Engine::ComplexV G1_hb, G2_hb;
    auto hb_prev = Engine::load(&Hb_lm1[(n - 1) * lanes]);

    if ((l - 1) == pl) {
      G1_hb = hb_prev;
      G2_hb = hb_prev;
    } else {
      auto term_hb = m_lm1 * hb_prev;
      G1_hb = term_hb - (m_l * d1_m1_n);
      auto d3_m1_n = Engine::load(&D3_mlxlM1[n * lanes]);
      G2_hb = term_hb - (m_l * d3_m1_n);
    }

    auto Temp_hb = Q_curr * G1_hb;
    auto Num_hb = (G2_hb * d1_n_curr) - (Temp_hb * d3_n);
    auto Denom_hb = G2_hb - Temp_hb;
    auto Hb_curr = Num_hb / Denom_hb;
    Engine::store(Hb_curr, &Hb_l[(n - 1) * lanes]);
  }
}

template <typename FloatType, typename Engine, typename XGetter, typename MGetter>
void calcScattCoeffsKernel(
    int nmax,
    int L,
    int pl,
    XGetter get_x,
    MGetter get_m,
    MieBuffers<FloatType, Engine>& buffers
) {
    int fl = (pl > 0) ? pl : 0;
    int lanes = Engine::Lanes();
    size_t stride = (nmax + 1) * lanes;

    // Aliases for buffers
    auto& D1_mlxl = buffers.D1;
    auto& D3_mlxl = buffers.D3;
    auto& D1_mlxlM1 = buffers.D1_prev;
    auto& D3_mlxlM1 = buffers.D3_prev;
    auto& PsiXL = buffers.Psi;
    auto& ZetaXL = buffers.Zeta;
    auto& Q = buffers.Q;
    auto& Ha = buffers.Ha;
    auto& Hb = buffers.Hb;
    auto& an = buffers.an;
    auto& bn = buffers.bn;

    // Layer fl (first layer)
    if (fl == pl) {
        // PEC
        auto zero = Engine::set(0.0);
        auto one = Engine::set(1.0);
        auto minus_one = Engine::set(-1.0);
        auto d1_val = Engine::make_complex(zero, minus_one);
        auto d3_val = Engine::make_complex(zero, one);
        
        for (int n = 0; n <= nmax; ++n) {
             Engine::store(d1_val, &D1_mlxl[n*lanes]);
             Engine::store(d3_val, &D3_mlxl[n*lanes]);
        }
    } else {
        auto x_fl = get_x(fl);
        auto m_fl = get_m(fl);
        auto zero = Engine::set(0.0);
        auto z1 = Engine::make_complex(x_fl, zero) * m_fl;

        evalDownwardD1<FloatType, Engine>(z1, D1_mlxl);
        evalUpwardD3<FloatType, Engine>(z1, D1_mlxl, D3_mlxl, PsiXL); // PsiXL used as temp buffer
    }
    
    // Ha, Hb for first layer
    for (int n = 0; n < nmax; ++n) {
        auto d1_np1 = Engine::load(&D1_mlxl[(n+1)*lanes]);
        Engine::store(d1_np1, &Ha[fl * stride + n * lanes]);
        Engine::store(d1_np1, &Hb[fl * stride + n * lanes]);
    }
    
    // Loop layers
    for (int l = fl + 1; l < L; ++l) {
        auto x_l = get_x(l);
        auto m_l = get_m(l);
        auto x_lm1 = get_x(l-1);
        auto m_lm1 = get_m(l-1);
        
        auto zero = Engine::set(0.0);
        auto z1 = Engine::make_complex(x_l, zero) * m_l;
        auto z2 = Engine::make_complex(x_lm1, zero) * m_l;

        evalDownwardD1<FloatType, Engine>(z1, D1_mlxl);
        evalUpwardD3<FloatType, Engine>(z1, D1_mlxl, D3_mlxl, PsiXL); // PsiXL temp
        
        evalDownwardD1<FloatType, Engine>(z2, D1_mlxlM1);
        evalUpwardD3<FloatType, Engine>(z2, D1_mlxlM1, D3_mlxlM1, PsiXL); // PsiXL temp

        computeLayerCoeffsHelper<FloatType, Engine>(
            nmax, l, pl, x_l, x_lm1, m_l, m_lm1, z1, z2, D1_mlxl.data(),
            D3_mlxl.data(), D1_mlxlM1.data(), D3_mlxlM1.data(), &Q[l * stride],
            &Ha[l * stride], &Hb[l * stride], &Ha[(l - 1) * stride],
            &Hb[(l - 1) * stride]);
    }
    
    // PsiXL, ZetaXL for outer layer
    auto x_L = get_x(L-1);
    auto zero = Engine::set(0.0);
    auto z_L = Engine::make_complex(x_L, zero);
    
    evalDownwardD1<FloatType, Engine>(z_L, D1_mlxl); // D1 temp
    evalUpwardPsi<FloatType, Engine>(z_L, D1_mlxl, PsiXL);
    evalUpwardD3<FloatType, Engine>(z_L, D1_mlxl, D3_mlxl, ZetaXL); // ZetaXL used as PsiZeta temp
    
    for (int n = 0; n <= nmax; ++n) {
        auto psi = Engine::load(&PsiXL[n*lanes]);
        auto psi_zeta = Engine::load(&ZetaXL[n*lanes]);
        auto zeta = psi_zeta / psi;
        Engine::store(zeta, &ZetaXL[n*lanes]);
    }
    
    // an, bn calculation
    auto m_L = get_m(L-1);
    auto one = Engine::set(1.0);
    auto one_c = Engine::make_complex(one, zero);
    auto zero_c = Engine::make_complex(zero, zero);

    for (int n = 0; n < nmax; ++n) {
        typename Engine::ComplexV an_val, bn_val;
        auto n_val = Engine::set(static_cast<FloatType>(n + 1));

        auto ha = Engine::load(&Ha[(L - 1) * stride + n * lanes]);
        auto hb = Engine::load(&Hb[(L - 1) * stride + n * lanes]);
        auto psi_np1 = Engine::load(&PsiXL[(n+1)*lanes]);
        auto zeta_np1 = Engine::load(&ZetaXL[(n+1)*lanes]);
        auto psi_n = Engine::load(&PsiXL[n*lanes]);
        auto zeta_n = Engine::load(&ZetaXL[n*lanes]);

        if (pl < (L - 1)) {
             computeAnBnBatch<FloatType, Engine>(
                  n_val, x_L, ha, hb, m_L,
                  psi_np1, zeta_np1, psi_n, zeta_n,
                  an_val, bn_val
             );
        } else {
             an_val = calc_an<FloatType, Engine>(
                  n_val, x_L, zero_c, one_c,
                  psi_np1, zeta_np1, psi_n, zeta_n
             );
             bn_val = psi_np1 / zeta_np1;
        }
        
        Engine::store(an_val, &an[n*lanes]);
        Engine::store(bn_val, &bn[n*lanes]);
    }
}

// ********************************************************************** //
// Calculates S1 - equation (25a)                                         //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
std::complex<FloatType> calc_S1(
    int n,
    std::complex<FloatType> an,
    std::complex<FloatType> bn,
    FloatType Pi,
    FloatType Tau) {
  return FloatType(n + n + 1) * (Pi * an + Tau * bn) / FloatType(n * n + n);
}

// ********************************************************************** //
// Calculates S2 - equation (25b) (it's the same as (25a), just switches  //
// Pi and Tau)                                                            //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
std::complex<FloatType> calc_S2(
    int n,
    std::complex<FloatType> an,
    std::complex<FloatType> bn,
    FloatType Pi,
    FloatType Tau) {
  return calc_S1(n, an, bn, Tau, Pi);
}



template <typename FloatType, typename Engine = ScalarEngine<FloatType>>
void calcQParams(
    int n,
    const std::complex<FloatType>& an,
    const std::complex<FloatType>& bn,
    const std::complex<FloatType>& an_next,
    const std::complex<FloatType>& bn_next,
    FloatType& Qext,
    FloatType& Qsca,
    FloatType& Qpr,
    std::complex<FloatType>& Qbktmp,
    int mode_n,
    int mode_type
) {
    const int n1 = n + 1;
    if (mode_n == Modes::kAll) {
      // Equation (27)
      Qext += (n1 + n1 + 1.0) * (an.real() + bn.real());
      // Equation (28)
      Qsca += (n1 + n1 + 1.0) *
               (an.real() * an.real() + an.imag() * an.imag() +
                bn.real() * bn.real() + bn.imag() * bn.imag());
      // Equation (29)
      Qpr += ((n1 * (n1 + 2.0) / (n1 + 1.0)) *
                   ((an * std::conj(an_next) + bn * std::conj(bn_next))
                        .real()) +
               ((n1 + n1 + 1.0) / (n1 * (n1 + 1.0))) *
                   (an * std::conj(bn)).real());
      // Equation (33)
      Qbktmp += (FloatType)(n1 + n1 + 1.0) * (1.0 - 2.0 * (n1 % 2)) *
                (an - bn);
    } else if (n1 == mode_n) {
      if (mode_type == Modes::kElectric || mode_type == Modes::kAll) {
        Qext += (n1 + n1 + 1.0) * (an.real());
        Qsca += (n1 + n1 + 1.0) * (an.real() * an.real() +
                                    an.imag() * an.imag());
        Qpr += std::nan("");
        Qbktmp +=
            (FloatType)(n1 + n1 + 1.0) * (1.0 - 2.0 * (n1 % 2)) * (an);
      }
      if (mode_type == Modes::kMagnetic || mode_type == Modes::kAll) {
        Qext += (n1 + n1 + 1.0) * (bn.real());
        Qsca += (n1 + n1 + 1.0) * (bn.real() * bn.real() +
                                    bn.imag() * bn.imag());
        Qpr += std::nan("");
        Qbktmp +=
            (FloatType)(n1 + n1 + 1.0) * (1.0 - 2.0 * (n1 % 2)) * (bn);
      }
    }
}

template <typename FloatType, typename Engine = ScalarEngine<FloatType>>
void calcSParams(
    int n,
    const std::complex<FloatType>& an,
    const std::complex<FloatType>& bn,
    const FloatType& Pi,
    const FloatType& Tau,
    std::complex<FloatType>& S1,
    std::complex<FloatType>& S2,
    int mode_n,
    int mode_type
) {
    const int n1 = n + 1;
    if (mode_n == Modes::kAll) {
        S1 += calc_S1<FloatType>(n1, an, bn, Pi, Tau);
        S2 += calc_S2<FloatType>(n1, an, bn, Pi, Tau);
    } else if (n1 == mode_n) {
      if (mode_type == Modes::kElectric || mode_type == Modes::kAll) {
          S1 += calc_S1<FloatType>(n1, an, std::complex<FloatType>(0), Pi, Tau);
          S2 += calc_S2<FloatType>(n1, an, std::complex<FloatType>(0), Pi, Tau);
      }
      if (mode_type == Modes::kMagnetic || mode_type == Modes::kAll) {
          S1 += calc_S1<FloatType>(n1, std::complex<FloatType>(0), bn, Pi, Tau);
          S2 += calc_S2<FloatType>(n1, std::complex<FloatType>(0), bn, Pi, Tau);
      }
    }
}

// class implementation

// ********************************************************************** //
// Returns previously calculated Qext                                     //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
template <typename outputType>
outputType MultiLayerMie<FloatType, Engine>::GetQext() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return static_cast<outputType>(Qext_);
}

// ********************************************************************** //
// Returns previously calculated Qabs                                     //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
template <typename outputType>
outputType MultiLayerMie<FloatType, Engine>::GetQabs() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return static_cast<outputType>(Qabs_);
}

// ********************************************************************** //
// Returns previously calculated Qsca                                     //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
template <typename outputType>
outputType MultiLayerMie<FloatType, Engine>::GetQsca() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return static_cast<outputType>(Qsca_);
}

// ********************************************************************** //
// Returns previously calculated Qbk                                      //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
template <typename outputType>
outputType MultiLayerMie<FloatType, Engine>::GetQbk() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return static_cast<outputType>(Qbk_);
}

// ********************************************************************** //
// Returns previously calculated Qpr                                      //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
template <typename outputType>
outputType MultiLayerMie<FloatType, Engine>::GetQpr() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return static_cast<outputType>(Qpr_);
}

// ********************************************************************** //
// Returns previously calculated asymmetry factor                         //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
template <typename outputType>
outputType MultiLayerMie<FloatType, Engine>::GetAsymmetryFactor() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return static_cast<outputType>(asymmetry_factor_);
}

// ********************************************************************** //
// Returns previously calculated Albedo                                   //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
template <typename outputType>
outputType MultiLayerMie<FloatType, Engine>::GetAlbedo() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return static_cast<outputType>(albedo_);
}

// ********************************************************************** //
// Returns previously calculated S1                                       //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
std::vector<std::complex<FloatType>> MultiLayerMie<FloatType, Engine>::GetS1() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return S1_;
}

// ********************************************************************** //
// Returns previously calculated S2                                       //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
std::vector<std::complex<FloatType>> MultiLayerMie<FloatType, Engine>::GetS2() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return S2_;
}

// ********************************************************************** //
// Modify scattering (theta) angles                                       //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
void MultiLayerMie<FloatType, Engine>::SetAngles(const std::vector<FloatType>& angles) {
  MarkUncalculated();
  theta_ = angles;
}

// ********************************************************************** //
// Modify size of all layers                                             //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
void MultiLayerMie<FloatType, Engine>::SetLayersSize(
    const std::vector<FloatType>& layer_size) {
  MarkUncalculated();
  size_param_.clear();
  FloatType prev_layer_size = 0.0;
  for (auto curr_layer_size : layer_size) {
    if (curr_layer_size <= 0.0)
      throw std::invalid_argument("Size parameter should be positive!");
    if (prev_layer_size > curr_layer_size)
      throw std::invalid_argument(
          "Size parameter for next layer should be larger than the previous "
          "one!");
    prev_layer_size = curr_layer_size;
    size_param_.push_back(curr_layer_size);
  }
}

// ********************************************************************** //
// Modify refractive index of all layers                                  //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
void MultiLayerMie<FloatType, Engine>::SetLayersIndex(
    const std::vector<std::complex<FloatType>>& index) {
  MarkUncalculated();
  refractive_index_ = index;
}

// ********************************************************************** //


// ********************************************************************** //
// Modify index of PEC layer                                              //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
void MultiLayerMie<FloatType, Engine>::SetPECLayer(int layer_position) {
  MarkUncalculated();
  if (layer_position < 0 && layer_position != -1)
    throw std::invalid_argument("Error! Layers are numbered from 0!");
  PEC_layer_position_ = layer_position;
}

// ********************************************************************** //
// Set maximun number of terms to be used                                 //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
void MultiLayerMie<FloatType, Engine>::SetMaxTerms(int nmax) const {
  MarkUncalculated();
  nmax_preset_ = nmax;
}

// ********************************************************************** //
// Get total size parameter of particle                                   //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
FloatType MultiLayerMie<FloatType, Engine>::GetSizeParameter() {
  if (size_param_.size() > 0)
    return size_param_.back();
  else
    return 0;
}

// ********************************************************************** //
// Mark uncalculated                                                      //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
void MultiLayerMie<FloatType, Engine>::MarkUncalculated() const {
  isExpCoeffsCalc_ = false;
  isScaCoeffsCalc_ = false;

  isMieCalculated_ = false;
}
// ********************************************************************** //
// Clear layer information                                                //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
void MultiLayerMie<FloatType, Engine>::ClearLayers() {
  MarkUncalculated();
  size_param_.clear();
  refractive_index_.clear();
}

// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
//                         Computational core
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //

template <typename FloatType>
unsigned int LeRu_near_field_cutoff(const std::complex<FloatType> zz) {
  std::complex<double> z = ConvertComplex<double>(zz);
  auto x = std::abs(z);
  return std::round(x + 11 * std::pow(x, (1.0 / 3.0)) + 1);
  //    return 10000;
}

// ********************************************************************** //
// Calculate calcNstop - equation (17)                                    //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
unsigned int MultiLayerMie<FloatType, Engine>::calcNstop(FloatType xL) const {
  unsigned int nmax = 0;
  // Wiscombe
  if (xL < size_param_.back())
    xL = size_param_.back();
  if (xL <= 8) {
    nmax = newround(xL + 4.0 * pow(xL, 1.0 / 3.0) + 1);
  } else if (xL <= 4200) {
    nmax = newround(xL + 4.05 * pow(xL, 1.0 / 3.0) + 2);
  } else {
    nmax = newround(xL + 4.0 * pow(xL, 1.0 / 3.0) + 2);
  }
  // Use Le Ru cutoff for near field, as a universal one.
  auto Nstop = nmie::LeRu_near_field_cutoff(std::complex<FloatType>(xL, 0)) + 1;
  if (Nstop > nmax)
    nmax = Nstop;
  return nmax;
}

// ********************************************************************** //
// Maximum number of terms required for the calculation                   //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
unsigned int MultiLayerMie<FloatType, Engine>::calcNmax(FloatType xL) const {
  const int pl = PEC_layer_position_;
  const unsigned int first_layer = (pl > 0) ? pl : 0;
  unsigned int ri, riM1, nmax = 0;
  const std::vector<FloatType>& x = size_param_;
  const std::vector<std::complex<FloatType>>& m = refractive_index_;
  nmax = calcNstop(xL);
  for (unsigned int i = first_layer; i < x.size(); i++) {
    if (static_cast<int>(i) >
        PEC_layer_position_)  // static_cast used to avoid warning
      ri = newround(cabs(x[i] * m[i]));
    else
      ri = 0;
    nmax = std::max(nmax, ri);
    // first layer is pec, if pec is present
    if ((i > first_layer) && (static_cast<int>(i - 1) > PEC_layer_position_))
      riM1 = newround(cabs(x[i - 1] * m[i]));
    else
      riM1 = 0;
    nmax = std::max(nmax, riM1);
  }
  nmax += 15;  // Final nmax value
#ifdef MULTI_PRECISION
  nmax += MULTI_PRECISION;  // TODO we may need to use more terms that this for
                            // MP computations.
#endif
  // nmax *= nmax;
  // printf("using nmax %i\n", nmax);
  return nmax;
}

// ********************************************************************** //
// Calculate an - equation (5)                                            //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
std::complex<FloatType> MultiLayerMie<FloatType, Engine>::calc_an(
    int n,
    FloatType XL,
    std::complex<FloatType> Ha,
    std::complex<FloatType> mL,
    std::complex<FloatType> PsiXL,
    std::complex<FloatType> ZetaXL,
    std::complex<FloatType> PsiXLM1,
    std::complex<FloatType> ZetaXLM1) const {
  return nmie::calc_an<FloatType, ScalarEngine<FloatType>, std::complex<FloatType>>(
      n, XL, Ha, mL, PsiXL, ZetaXL, PsiXLM1, ZetaXLM1);
}

// ********************************************************************** //
// Calculate bn - equation (6)                                            //
// ********************************************************************** //
template <typename FloatType, MathEngine Engine>
std::complex<FloatType> MultiLayerMie<FloatType, Engine>::calc_bn(
    int n,
    FloatType XL,
    std::complex<FloatType> Hb,
    std::complex<FloatType> mL,
    std::complex<FloatType> PsiXL,
    std::complex<FloatType> ZetaXL,
    std::complex<FloatType> PsiXLM1,
    std::complex<FloatType> ZetaXLM1) const {
  return nmie::calc_bn<FloatType, ScalarEngine<FloatType>, std::complex<FloatType>>(
      n, XL, Hb, mL, PsiXL, ZetaXL, PsiXLM1, ZetaXLM1);
}



//****************************************************************************
// This function calculates the logarithmic derivatives of the Riccati-Bessel
// functions (D1 and D3) for a complex argument (z).
// Equations (16a), (16b) and (18a) - (18d)
//
// Input parameters:
//   z: Complex argument to evaluate D1 and D3
//   nmax_: Maximum number of terms to calculate D1 and D3
//
// Output parameters:
//   D1, D3: Logarithmic derivatives of the Riccati-Bessel functions
//****************************************************************************
template <typename FloatType>
void calcD1D3(
    const std::complex<FloatType> z,
    int nmax,
    std::vector<std::complex<FloatType>>& D1,
    std::vector<std::complex<FloatType>>& D3) {
  std::vector<std::complex<FloatType>> PsiZeta(nmax + 1);
  evalDownwardD1<FloatType>(z, D1);
  evalUpwardD3<FloatType>(z, D1, D3, PsiZeta);
}

//*****************************************************************************
// This function calculates the Riccati-Bessel functions (Psi and Zeta) for a
// complex argument (z).
// Equations (20a) - (21b)
//
// Input parameters:
//   z: Complex argument to evaluate Psi and Zeta
//   nmax: Maximum number of terms to calculate Psi and Zeta
//
// Output parameters:
//   Psi, Zeta: Riccati-Bessel functions
//*****************************************************************************
template <typename FloatType>
void calcPsiZeta(
    std::complex<FloatType> z,
    int nmax,
    std::vector<std::complex<FloatType>>& Psi,
    std::vector<std::complex<FloatType>>& Zeta) {
  std::vector<std::complex<FloatType>> D1(nmax + 1), D3(nmax + 1),
      PsiZeta(nmax + 1);
  // First, calculate the logarithmic derivatives
  evalDownwardD1<FloatType>(z, D1);
  // Now, use the upward recurrence to calculate Psi equations (20ab)
  evalUpwardPsi<FloatType>(z, D1, Psi);
  // Now, use the upward recurrence to calculate Psi*Zeta equations (18ad)
  evalUpwardD3<FloatType>(z, D1, D3, PsiZeta);
  for (unsigned int i = 0; i < Zeta.size(); i++) {
    Zeta[i] = PsiZeta[i] / Psi[i];
  }
  //    evalUpwardZeta(z, D3, Zeta);
}

template <typename FloatType, MathEngine Engine>
void MultiLayerMie<FloatType, Engine>::calcPiTauAllTheta(
    const double from_Theta,
    const double to_Theta,
    std::vector<std::vector<FloatType>>& Pi,
    std::vector<std::vector<FloatType>>& Tau) const {
  const unsigned int perimeter_points = Pi.size();
  for (auto& val : Pi)
    val.resize(available_maximal_nmax_, static_cast<FloatType>(0.0));
  for (auto& val : Tau)
    val.resize(available_maximal_nmax_, static_cast<FloatType>(0.0));
  double delta_Theta =
      eval_delta<double>(perimeter_points, from_Theta, to_Theta);
  for (unsigned int i = 0; i < perimeter_points; i++) {
    auto Theta = static_cast<FloatType>(from_Theta + i * delta_Theta);
    // Calculate angular functions Pi and Tau
    calcPiTau(nmm::cos(Theta), Pi[i], Tau[i]);
  }
}

//*******************************************************************************
// This function calculates Pi and Tau for a given value of cos(Theta).
// Equations (26a) - (26c)
//
// Input parameters:
//   nmax_: Maximum number of terms to calculate Pi and Tau
//   nTheta: Number of scattering angles
//   Theta: Array containing all the scattering angles where the scattering
//          amplitudes will be calculated
//
// Output parameters:
//   Pi, Tau: Angular functions Pi and Tau, as defined in equations (26a) -
//   (26c)
//*******************************************************************************
template <typename FloatType, MathEngine Engine>
void MultiLayerMie<FloatType, Engine>::calcPiTau(const FloatType& costheta,
                                         std::vector<FloatType>& Pi,
                                         std::vector<FloatType>& Tau) const {
  int nmax = Pi.size();
  if (Pi.size() != Tau.size())
    throw std::invalid_argument(
        "Error! Pi and Tau vectors should have the same size!");

  //****************************************************//
  // Equations (26a) - (26c)                            //
  //****************************************************//
  // Initialize Pi and Tau
  Pi[0] = 1.0;  // n=1
  Tau[0] = costheta;
  // Calculate the actual values
  if (nmax > 1) {
    Pi[1] = 3 * costheta * Pi[0];  // n=2
    Tau[1] = 2 * costheta * Pi[1] - 3 * Pi[0];
    for (int i = 2; i < nmax; i++) {  // n=[3..nmax_]
      Pi[i] = ((i + i + 1) * costheta * Pi[i - 1] - (i + 1) * Pi[i - 2]) / i;
      Tau[i] = (i + 1) * costheta * Pi[i] - (i + 2) * Pi[i - 1];
    }
  }
}  // end of MultiLayerMie::calcPiTau(...)

//*****************************************************************************
// This function calculates vector spherical harmonics (eq. 4.50, p. 95 BH),
// required to calculate the near-field parameters.
//
// Input parameters:
//   Rho: Radial distance
//   Phi: Azimuthal angle
//   Theta: Polar angle
//   rn: Either the spherical Ricatti-Bessel function of first or third kind
//   Dn: Logarithmic derivative of rn
//   Pi, Tau: Angular functions Pi and Tau
//   n: Order of vector spherical harmonics
//
// Output parameters:
//   Mo1n, Me1n, No1n, Ne1n: Complex vector spherical harmonics
//*****************************************************************************
template <typename FloatType, MathEngine Engine>
template <typename evalType>
void MultiLayerMie<FloatType, Engine>::calcSpherHarm(
    const std::complex<evalType> Rho,
    const evalType Theta,
    const evalType Phi,
    const std::complex<evalType>& rn,
    const std::complex<evalType>& Dn,
    const evalType& Pi,
    const evalType& Tau,
    const evalType& n,
    std::vector<std::complex<evalType>>& Mo1n,
    std::vector<std::complex<evalType>>& Me1n,
    std::vector<std::complex<evalType>>& No1n,
    std::vector<std::complex<evalType>>& Ne1n) const {
  // using eq 4.50 in BH
  std::complex<evalType> c_zero(0.0, 0.0);

  //    using nmm::sin;
  //    using nmm::cos;
  auto sin_Phi = sin_t(Phi);
  auto cos_Phi = cos_t(Phi);
  auto sin_Theta = sin(Theta);
  Mo1n[0] = c_zero;
  Mo1n[1] = cos_Phi * Pi * rn / Rho;
  Mo1n[2] = -sin_Phi * Tau * rn / Rho;

  Me1n[0] = c_zero;
  Me1n[1] = -sin_Phi * Pi * rn / Rho;
  Me1n[2] = -cos_Phi * Tau * rn / Rho;

  No1n[0] = sin_Phi * (n * n + n) * sin_Theta * Pi * rn / Rho / Rho;
  No1n[1] = sin_Phi * Tau * Dn * rn / Rho;
  No1n[2] = cos_Phi * Pi * Dn * rn / Rho;

  Ne1n[0] = cos_Phi * (n * n + n) * sin_Theta * Pi * rn / Rho / Rho;
  Ne1n[1] = cos_Phi * Tau * Dn * rn / Rho;
  Ne1n[2] = -sin_Phi * Pi * Dn * rn / Rho;
}  // end of MultiLayerMie::calcSpherHarm(...)

//********************************************************************************
// This function calculates the scattering coefficients required to calculate
// both the near- and far-field parameters.
//
// Input parameters:
//   L: Number of layers
//   pl: Index of PEC layer. If there is none just send -1
//   x: Array containing the size parameters of the layers [0..L-1]
//   m: Array containing the relative refractive indexes of the layers [0..L-1]
//   nmax: Maximum number of multipolar expansion terms to be used for the
//         calculations. Only use it if you know what you are doing, otherwise
//         set this parameter to -1 and the function will calculate it.
//
// Output parameters:
//   an, bn: Complex scattering amplitudes
//
// Return value:
//   Number of multipolar expansion terms used for the calculations
//********************************************************************************
template <typename FloatType, MathEngine Engine>
void MultiLayerMie<FloatType, Engine>::calcScattCoeffs() const {
  isScaCoeffsCalc_ = false;
  an_.clear();
  bn_.clear();

  const std::vector<FloatType>& x = size_param_;
  const std::vector<std::complex<FloatType>>& m = refractive_index_;
  const int& pl = PEC_layer_position_;
  const int L = refractive_index_.size();

  //************************************************************************//
  // Calculate the index of the first layer. It can be either 0 (default)   //
  // or the index of the outermost PEC layer. In the latter case all layers //
  // below the PEC are discarded.                                           //
  // ***********************************************************************//
  if (nmax_preset_ <= 0)
    nmax_ = calcNmax();
  else
    nmax_ = nmax_preset_;

  //**************************************************************************//
  // Note that since Fri, Nov 14, 2014 all arrays start from 0 (zero), which  //
  // means that index = layer number - 1 or index = n - 1. The only exception //
  // are the arrays for representing D1, D3 and Q because they need a value   //
  // for the index 0 (zero), hence it is important to consider this shift     //
  // between different arrays. The change was done to optimize memory usage.  //
  //**************************************************************************//
  
  size_t lanes = Engine::Lanes();
  MieBuffers<FloatType, Engine> buffers;
  buffers.resize(nmax_, L, 0);
  
  auto get_x = [&](int l) { return Engine::set(x[l]); };
  auto get_m = [&](int l) { 
    auto m_val = m[l];
    return Engine::make_complex(
        Engine::set(m_val.real()), 
        Engine::set(m_val.imag())
    ); 
  };
  
  calcScattCoeffsKernel<FloatType, Engine>(
    nmax_, L, pl, get_x, get_m,
    buffers
  );
  
  an_.resize(nmax_);
  bn_.resize(nmax_);
  for(int n=0; n<nmax_; ++n) {
    an_[n] = buffers.an[n * lanes];
    bn_[n] = buffers.bn[n * lanes];
  }

  FloatType a0 = 0, b0 = 0;

  for (int n = 0; n < nmax_; n++) {
    if (n == 0) {
      a0 = cabs(an_[0]);
      b0 = cabs(bn_[0]);
    }
    if (n == nmax_ - 1 && nmax_preset_ <= 0 &&
        (cabs(an_[n]) / a0 > convergence_threshold_ &&
         cabs(bn_[n]) / b0 > convergence_threshold_)) {
      std::cout << "Failed to converge in Mie series for nmax=" << nmax_
                << std::endl;
      std::cout << "convergence threshold: " << convergence_threshold_
                << std::endl;
      std::cout << "Mie series a[nmax]/a[1]:" << cabs(an_[n]) / a0
                << " and b[nmax]/b[1]:" << cabs(bn_[n]) / b0 << std::endl;
    }

    if (nmm::isnan(an_[n].real()) || nmm::isnan(an_[n].imag()) ||
        nmm::isnan(bn_[n].real()) || nmm::isnan(bn_[n].imag())) {
      std::cout
          << "nmax value was changed due to unexpected error!!! New values is "
          << n << " (was " << nmax_ << ")" << std::endl;
      nmax_ = n;
      an_.resize(nmax_);
      bn_.resize(nmax_);
      break;
    }

  }  // end of for an and bn terms
  isScaCoeffsCalc_ = true;
}  // end of MultiLayerMie::calcScattCoeffs()

template <typename FloatType, typename Engine>
void sumMieSeriesKernel(const int nmax,
                        const typename Engine::RealV vnmax,
                        const typename Engine::RealV x,
                        const std::complex<FloatType>* an_ptr,
                        const std::complex<FloatType>* bn_ptr,
                        const std::vector<FloatType>& theta,
                        MieBuffers<FloatType, Engine>& buffers,
                        std::vector<typename Engine::ComplexV>& S1,
                        std::vector<typename Engine::ComplexV>& S2) {
  using ComplexV = typename Engine::ComplexV;
  using RealV = typename Engine::RealV;

  auto zero = Engine::set(0.0);
  auto one = Engine::set(1.0);

  auto Qext_sum = zero;
  auto Qsca_sum = zero;
  auto Qpr_sum = zero;
  auto Qbk_sum = Engine::make_complex(zero, zero);

  size_t num_angles = theta.size();
  for (size_t k = 0; k < num_angles; ++k) {
    S1[k] = Engine::make_complex(zero, zero);
    S2[k] = Engine::make_complex(zero, zero);
  }

  std::vector<ComplexV> pi_prev(num_angles, Engine::make_complex(zero, zero));
  std::vector<ComplexV> pi_curr(num_angles, Engine::make_complex(one, zero));
  std::vector<RealV> mu(num_angles);
  for (size_t k = 0; k < num_angles; ++k) {
    mu[k] = Engine::set(nmm::cos(theta[k]));
  }

  ComplexV an_prev = Engine::make_complex(zero, zero);
  ComplexV bn_prev = Engine::make_complex(zero, zero);

  for (int n = 1; n < nmax; ++n) {
    auto vn = Engine::set(static_cast<FloatType>(n));
    auto mask = Engine::le(vn, vnmax);

    size_t offset = (n - 1) * Engine::Lanes();
    auto an = Engine::load(&an_ptr[offset]);
    auto bn = Engine::load(&bn_ptr[offset]);

    an = Engine::select(mask, an, Engine::make_complex(zero, zero));
    bn = Engine::select(mask, bn, Engine::make_complex(zero, zero));

    auto n2p1 = (vn + vn) + one;
    auto n_plus_1 = vn + one;
    auto n_sq_plus_n = vn * n_plus_1;

    auto an_plus_bn = an + bn;
    auto term_ext = n2p1 * Engine::get_real(an_plus_bn);
    Qext_sum = Qext_sum + term_ext;

    auto an_sq = (Engine::get_real(an) * Engine::get_real(an)) +
                 (Engine::get_imag(an) * Engine::get_imag(an));
    auto bn_sq = (Engine::get_real(bn) * Engine::get_real(bn)) +
                 (Engine::get_imag(bn) * Engine::get_imag(bn));
    auto term_sca = n2p1 * (an_sq + bn_sq);
    Qsca_sum = Qsca_sum + term_sca;

    RealV sign = (n % 2 == 0) ? one : Engine::set(-1.0);
    auto an_minus_bn = an - bn;
    auto term_bk = Engine::make_complex(n2p1 * sign, zero) * an_minus_bn;
    Qbk_sum = Qbk_sum + term_bk;

    if (n > 1) {
      auto nm1 = vn - one;
      auto factor1 = (nm1 * n_plus_1) / vn;

      auto re_aa = (Engine::get_real(an_prev) * Engine::get_real(an)) +
                   (Engine::get_imag(an_prev) * Engine::get_imag(an));
      auto re_bb = (Engine::get_real(bn_prev) * Engine::get_real(bn)) +
                   (Engine::get_imag(bn_prev) * Engine::get_imag(bn));

      auto term_pr1 = factor1 * (re_aa + re_bb);
      Qpr_sum = Qpr_sum + term_pr1;
    }

    auto factor2 = n2p1 / n_sq_plus_n;
    auto re_ab = (Engine::get_real(an) * Engine::get_real(bn)) +
                 (Engine::get_imag(an) * Engine::get_imag(bn));
    auto term_pr2 = factor2 * re_ab;
    Qpr_sum = Qpr_sum + term_pr2;

    an_prev = an;
    bn_prev = bn;

    if (num_angles > 0) {
      auto factor = factor2;

      for (size_t k = 0; k < num_angles; ++k) {
        ComplexV pi_next = Engine::make_complex(zero, zero);
        ComplexV tau_n = Engine::make_complex(zero, zero);

        if (n == 1) {
          pi_next = Engine::make_complex(one, zero);
          tau_n = Engine::make_complex(mu[k], zero);
        } else {
          auto nm1 = vn - one;
          auto n2m1 = (vn + vn) - one;

          auto t1 = n2m1 * (mu[k] * Engine::get_real(pi_curr[k]));
          auto t2 = vn * Engine::get_real(pi_prev[k]);
          auto pi_val = (t1 - t2) / nm1;
          pi_next = Engine::make_complex(pi_val, zero);

          auto t3 = vn * (mu[k] * pi_val);
          auto t4 = n_plus_1 * Engine::get_real(pi_curr[k]);
          auto tau_val = t3 - t4;
          tau_n = Engine::make_complex(tau_val, zero);

          pi_prev[k] = pi_curr[k];
          pi_curr[k] = pi_next;
        }

        auto term_S1 = (an * pi_curr[k]) + (bn * tau_n);
        S1[k] = S1[k] + (Engine::make_complex(factor, zero) * term_S1);

        auto term_S2 = (an * tau_n) + (bn * pi_curr[k]);
        S2[k] = S2[k] + (Engine::make_complex(factor, zero) * term_S2);
      }
    }
  }

  // Finalize results
  auto two = Engine::set(2.0);
  auto four = Engine::set(4.0);
  auto x2 = x * x;
  auto norm = two / x2;

  buffers.Qext = Qext_sum * norm;
  buffers.Qsca = Qsca_sum * norm;
  buffers.Qpr = buffers.Qext - (four * Qpr_sum / x2);
  buffers.Qabs = buffers.Qext - buffers.Qsca;

  auto qbk_mag_sq = (Engine::get_real(Qbk_sum) * Engine::get_real(Qbk_sum)) +
                    (Engine::get_imag(Qbk_sum) * Engine::get_imag(Qbk_sum));
  buffers.Qbk = qbk_mag_sq / x2;

  auto threshold = Engine::set(1e-12);
  auto mask_sca = Engine::gt(buffers.Qsca, threshold);
  auto mask_ext = Engine::gt(buffers.Qext, threshold);

  buffers.g = Engine::select(mask_sca,
                             (buffers.Qext - buffers.Qpr) / buffers.Qsca, zero);
  buffers.albedo = Engine::select(mask_ext, buffers.Qsca / buffers.Qext, zero);

  // Store results to _res fields
  buffers.Qext_res = Engine::to_scalar(buffers.Qext);
  buffers.Qsca_res = Engine::to_scalar(buffers.Qsca);
  buffers.Qabs_res = Engine::to_scalar(buffers.Qabs);
  buffers.Qbk_res = Engine::to_scalar(buffers.Qbk);
  buffers.Qpr_res = Engine::to_scalar(buffers.Qpr);
  buffers.g_res = Engine::to_scalar(buffers.g);
  buffers.albedo_res = Engine::to_scalar(buffers.albedo);
}

template <typename FloatType, typename Engine>
void finalizeMieResults(const typename Engine::RealV x,
                        const typename Engine::RealV Qext_sum,
                        const typename Engine::RealV Qsca_sum,
                        const typename Engine::RealV Qpr_sum,
                        const typename Engine::ComplexV Qbk_sum,
                        typename Engine::RealV& Qext,
                        typename Engine::RealV& Qsca,
                        typename Engine::RealV& Qabs,
                        typename Engine::RealV& Qbk,
                        typename Engine::RealV& Qpr,
                        typename Engine::RealV& g,
                        typename Engine::RealV& Albedo) {
  auto two = Engine::set(2.0);
  auto four = Engine::set(4.0);
  auto x2 = x * x;
  auto norm = two / x2;

  Qext = Qext_sum * norm;
  Qsca = Qsca_sum * norm;
  Qpr = Qext - (four * Qpr_sum / x2);
  Qabs = Qext - Qsca;

  auto qbk_mag_sq = (Engine::get_real(Qbk_sum) * Engine::get_real(Qbk_sum)) +
                    (Engine::get_imag(Qbk_sum) * Engine::get_imag(Qbk_sum));
  Qbk = qbk_mag_sq / x2;

  auto zero = Engine::set(0.0);
  auto threshold = Engine::set(1e-12);

  auto mask_sca = Engine::gt(Qsca, threshold);
  auto mask_ext = Engine::gt(Qext, threshold);

  g = Engine::select(mask_sca, (Qext - Qpr) / Qsca, zero);
  Albedo = Engine::select(mask_ext, Qsca / Qext, zero);
}

//*******************************************************************************
// This function calculates the actual scattering parameters and amplitudes
//
// Input parameters:
//   L: Number of layers
//   pl: Index of PEC layer. If there is none just send -1
//   x: Array containing the size parameters of the layers [0..L-1]
//   m: Array containing the relative refractive indexes of the layers [0..L-1]
//   nTheta: Number of scattering angles
//   Theta: Array containing all the scattering angles where the scattering
//          amplitudes will be calculated
//   nmax_: Maximum number of multipolar expansion terms to be used for the
//         calculations. Only use it if you know what you are doing, otherwise
//         set this parameter to -1 and the function will calculate it
//
// Output parameters:
//   Qext: Efficiency factor for extinction
//   Qsca: Efficiency factor for scattering
//   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)
//   Qbk: Efficiency factor for backscattering
//   Qpr: Efficiency factor for the radiation pressure
//   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)
//   Albedo: Single scattering albedo (Albedo = Qsca/Qext)
//   S1, S2: Complex scattering amplitudes
//
// Return value:
//   Number of multipolar expansion terms used for the calculations
//*******************************************************************************
template <typename FloatType, MathEngine Engine>
void MultiLayerMie<FloatType, Engine>::RunMieCalculationStateless(
    MieBuffers<FloatType, Engine>& buffers) const {
  if (size_param_.size() != refractive_index_.size())
    throw std::invalid_argument(
        "Each size parameter should have only one index!");
  if (size_param_.size() == 0)
    throw std::invalid_argument("Initialize model first!");

  const std::vector<FloatType>& x = size_param_;

  // Initialize the scattering parameters
  buffers.Qext_res = 0.0;
  buffers.Qsca_res = 0.0;
  buffers.Qabs_res = 0.0;
  buffers.Qbk_res = 0.0;
  buffers.Qpr_res = 0.0;

  buffers.g_res = 0.0;
  buffers.albedo_res = 0.0;

  // Initialize the scattering amplitudes
  // buffers.S1_res and S2_res are scalar vectors.
  // We need SIMD vectors for the kernel.
  std::vector<typename Engine::ComplexV> S1_simd(theta_.size());
  std::vector<typename Engine::ComplexV> S2_simd(theta_.size());

  // Prepare an and bn for SIMD kernel
  // sumMieSeriesKernel expects interleaved/SIMD-ready data
  // We need to broadcast an_ and bn_ to all lanes
  size_t lanes = Engine::Lanes();
  std::vector<std::complex<FloatType>> an_simd(nmax_ * lanes);
  std::vector<std::complex<FloatType>> bn_simd(nmax_ * lanes);
  
  for (int n = 0; n < nmax_; ++n) {
      for (size_t l = 0; l < lanes; ++l) {
          an_simd[n * lanes + l] = an_[n];
          bn_simd[n * lanes + l] = bn_[n];
      }
  }

  sumMieSeriesKernel<FloatType, Engine>(
      nmax_, 
      Engine::set(static_cast<FloatType>(nmax_)), 
      Engine::set(x.back()), 
      an_simd.data(), bn_simd.data(),
      theta_, buffers, S1_simd, S2_simd);

  // Copy back S1_simd to buffers.S1_res
  buffers.S1_res.resize(theta_.size());
  buffers.S2_res.resize(theta_.size());
  for (size_t i = 0; i < theta_.size(); ++i) {
    buffers.S1_res[i] = std::complex<FloatType>(
        Engine::to_scalar(Engine::get_real(S1_simd[i])), 
        Engine::to_scalar(Engine::get_imag(S1_simd[i])));
    buffers.S2_res[i] = std::complex<FloatType>(
        Engine::to_scalar(Engine::get_real(S2_simd[i])), 
        Engine::to_scalar(Engine::get_imag(S2_simd[i])));
  }
}

template <typename FloatType, MathEngine Engine>
void MultiLayerMie<FloatType, Engine>::RunMieCalculation() const {
  if (size_param_.size() != refractive_index_.size())
    throw std::invalid_argument(
        "Each size parameter should have only one index!");
  if (size_param_.size() == 0)
    throw std::invalid_argument("Initialize model first!");

  // MarkUncalculated();

  // Calculate scattering coefficients
  if (!isScaCoeffsCalc_)
    calcScattCoeffs();

  MieBuffers<FloatType, Engine> buffers;
  buffers.updateSize(nmax_, size_param_.size(), theta_.size());

  RunMieCalculationStateless(buffers);

  Qext_ = buffers.Qext_res;
  Qsca_ = buffers.Qsca_res;
  Qabs_ = buffers.Qabs_res;
  Qbk_ = buffers.Qbk_res;
  Qpr_ = buffers.Qpr_res;
  asymmetry_factor_ = buffers.g_res;
  albedo_ = buffers.albedo_res;
  S1_ = buffers.S1_res;
  S2_ = buffers.S2_res;

  isMieCalculated_ = true;
}

}  // end of namespace nmie
#endif  // SRC_NMIE_BASIC_HPP_
