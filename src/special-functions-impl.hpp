#ifndef SRC_SPECIAL_FUNCTIONS_IMPL_HPP_
#define SRC_SPECIAL_FUNCTIONS_IMPL_HPP_
//******************************************************************************
//    Copyright (C) 2009-2022  Ovidio Pena <ovidio@bytesfall.com>
//    Copyright (C) 2013-2022  Konstantin Ladutenko <kostyfisik@gmail.com>
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
//******************************************************************************
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "nmie-precision.hpp"
#include "nmie.hpp"

namespace nmie {

//******************************************************************************
// Note, that Kapteyn seems to be too optimistic (at least by 3 digits
// in some cases) for forward recurrence, see D1test with WYang_data
//******************************************************************************
template <typename FloatType, typename Engine = ScalarEngine<FloatType>, typename ComplexType = std::complex<FloatType>>
typename Engine::RealV evalKapteynNumberOfLostSignificantDigits(const int ni, const ComplexType zz) {
  auto z = zz;
  auto n_val = Engine::set(static_cast<FloatType>(ni));
  auto one_val = Engine::set(1.0);
  auto zero_val = Engine::set(0.0);
  auto one = Engine::make_complex(one_val, zero_val);
  
  auto n_complex = Engine::make_complex(n_val, zero_val);
  auto z_div_n = z / n_complex;
  
  auto z_div_n_sq = z_div_n * z_div_n;
  auto one_minus_sq = one - z_div_n_sq;
  auto sqrt_term = Engine::sqrt(one_minus_sq);
  
  auto log_z_div_n = Engine::log(z_div_n);
  auto one_plus_sqrt = one + sqrt_term;
  auto log_one_plus_sqrt = Engine::log(one_plus_sqrt);
  
  auto term = (log_z_div_n + sqrt_term) - log_one_plus_sqrt;
  auto real_term = Engine::get_real(term);
  auto n_times_real = n_val * real_term;
  
  auto imag_z = Engine::get_imag(z);
  auto abs_imag_z = Engine::abs(imag_z);
  auto log_2 = Engine::set(std::log(2.0));
  
  auto num = (abs_imag_z - log_2) - n_times_real;
  auto log_10 = Engine::set(std::log(10.0));
  
  auto res = num / log_10;
  auto half = Engine::set(0.5);
  auto res_plus_half = res + half;
  auto rounded = Engine::floor(res_plus_half);
  
  return rounded;
}

//******************************************************************************
template <typename FloatType, typename Engine = ScalarEngine<FloatType>, typename ComplexType = std::complex<FloatType>>
int getNStar(int nmax, ComplexType z, const int valid_digits) {
  if (nmax == 0) nmax = 1;
  int nstar = nmax;
  (void)valid_digits;
  
  auto forwardLoss = evalKapteynNumberOfLostSignificantDigits<FloatType, Engine, ComplexType>(nmax, z);
  
  auto re = Engine::get_real(z);
  auto im = Engine::get_imag(z);
  auto abs_z_val = Engine::sqrt((re * re) + (im * im));
  
  auto pow_third = Engine::pow(abs_z_val, Engine::set(1.0/3.0));
  auto four_pow = Engine::set(15.0) * pow_third;
  auto max_val = Engine::max(four_pow, Engine::set(5.0));
  auto increment_v = Engine::ceil(max_val);
  
  int increment = static_cast<int>(Engine::reduce_max(increment_v));
  
  auto max_abs_z_v = Engine::ceil(abs_z_val);
  int max_abs_z = static_cast<int>(Engine::reduce_max(max_abs_z_v));

  if (nstar < max_abs_z) nstar = max_abs_z;
  nstar += increment;

  int iter = 0;
  while (iter++ < 100) {
    auto backwardLoss = evalKapteynNumberOfLostSignificantDigits<FloatType, Engine, ComplexType>(nstar, z);
    auto diff = backwardLoss - forwardLoss;
    
    auto max_diff = Engine::reduce_max(diff);
    if (max_diff > valid_digits) break;
    
    nstar += increment;
  }

  return nstar;
}

//******************************************************************************
// Custom implementation of complex cot function to avoid overflow
// if Im(z) < 0, then it evaluates cot(z) as conj(cot(conj(z)))
// see Eqs. 10-12 of [1] for details.
// [1]H. Du, Mie-Scattering Calculation, Appl. Opt. 43, 1951 (2004).
//******************************************************************************
template <typename FloatType, typename Engine = ScalarEngine<FloatType>, typename ComplexType>
ComplexType complex_cot(const ComplexType z) {
  auto Remx = Engine::get_real(z);
  auto Immx = Engine::get_imag(z);
  
  auto sgn = Engine::sign(Immx);
  
  auto minus_two = Engine::set(-2.0);
  auto arg = minus_two * sgn * Immx;
  auto exp_val = Engine::exp(arg);
  
  auto tan_val = Engine::tan(Remx);
  
  auto exp_tan = exp_val * tan_val;
  auto a = tan_val - exp_tan;
  
  auto one = Engine::set(1.0);
  auto b = one + exp_val;
  
  auto minus_one = Engine::set(-1.0);
  auto c = minus_one + exp_val;
  
  auto d = tan_val + exp_tan;
  
  auto c2 = c * c;
  auto d2 = d * d;
  auto denom = c2 + d2;
  
  auto ac = a * c;
  auto bd = b * d;
  auto num_re = ac + bd;
  auto re = num_re / denom;
  
  auto bc = b * c;
  auto ad = a * d;
  auto diff = bc - ad;
  auto num_im = sgn * diff;
  auto im = num_im / denom;
  
  return Engine::make_complex(re, im);
}

//******************************************************************************
// Forward iteration for evaluation of ratio of the Riccati–Bessel functions
//******************************************************************************
template <typename FloatType>
void evalForwardR(const std::complex<FloatType> z,
                  std::vector<std::complex<FloatType>>& r) {
  if (r.size() < 1)
    throw std::invalid_argument(
        "We need non-zero evaluations of ratio of the Riccati–Bessel "
        "functions.\n");
  // r0 = cot(z)
  //  r[0] = nmm::cos(z)/nmm::sin(z);
  r[0] = complex_cot<FloatType>(z);
  for (unsigned int n = 0; n < r.size() - 1; n++) {
    r[n + 1] = static_cast<FloatType>(1) /
               (                                          //
                   static_cast<FloatType>(2 * n + 1) / z  //
                   - r[n]                                 //
               );
  }
}

//******************************************************************************
// Backward iteration for evaluation of ratio of the Riccati–Bessel functions
//******************************************************************************
template <typename FloatType>
void evalBackwardR(const std::complex<FloatType> z,
                   std::vector<std::complex<FloatType>>& r) {
  if (r.size() < 1)
    throw std::invalid_argument(
        "We need non-zero evaluations of ratio of the Riccati-Bessel "
        "functions.\n");
  int nmax = r.size() - 1;
  auto tmp = static_cast<FloatType>(2 * nmax + 1) / z;
  r[nmax] = tmp;
  for (int n = nmax - 1; n >= 0; n--) {
    r[n] = -static_cast<FloatType>(1) / r[n + 1] +
           static_cast<FloatType>(2 * n + 1) / z;
  }
  r[0] = complex_cot<FloatType>(z);
}

//******************************************************************************
template <typename FloatType>
void convertRtoD1(const std::complex<FloatType> z,
                  std::vector<std::complex<FloatType>>& r,
                  std::vector<std::complex<FloatType>>& D1) {
  if (D1.size() > r.size())
    throw std::invalid_argument(
        "Not enough elements in array of ratio of the Riccati–Bessel functions "
        "to convert it into logarithmic derivative array.\n");
  std::complex<FloatType> D_old;
  for (unsigned int n = 0; n < D1.size(); n++) {
    D_old = D1[n];
    D1[n] = r[n] - static_cast<FloatType>(n) / z;
    //    std::cout << "D1old - D1[n] :" << D_old-D1[n] << '\n';
  }
}

//******************************************************************************
template <typename FloatType>
void evalForwardD(const std::complex<FloatType> z,
                  std::vector<std::complex<FloatType>>& D) {
  int nmax = D.size();
  FloatType one = static_cast<FloatType>(1);
  for (int n = 1; n < nmax; n++) {
    FloatType nf = static_cast<FloatType>(n);
    D[n] = one / (nf / z - D[n - 1]) - nf / z;
  }
}

//******************************************************************************
template <typename FloatType>
void evalForwardD1(const std::complex<FloatType> z,
                   std::vector<std::complex<FloatType>>& D) {
  if (D.size() < 1)
    throw std::invalid_argument("Should have a leas one element!\n");
  D[0] = nmm::cos(z) / nmm::sin(z);
  evalForwardD(z, D);
}

//  template <typename FloatType>
//  void MultiLayerMie<FloatType>::calcD1D3(const std::complex<FloatType> z,
//                               std::vector<std::complex<FloatType> > &D1,
//                               std::vector<std::complex<FloatType> > &D3) {
//
//    // if (cabs(D1[0]) > 1.0e15) {
//    //   throw std::invalid_argument("Unstable D1! Please, try to change input
//    parameters!\n");
//    // //printf("Warning: Potentially unstable D1! Please, try to change input
//    parameters!\n");
//    // }
//
//    // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
//    PsiZeta_[0] = static_cast<FloatType>(0.5)*(static_cast<FloatType>(1.0) -
//    std::complex<FloatType>(nmm::cos(2.0*z.real()), nmm::sin(2.0*z.real()))
//                 *static_cast<FloatType>(nmm::exp(-2.0*z.imag())));
//    D3[0] = std::complex<FloatType>(0.0, 1.0);
//
//    for (int n = 1; n <= nmax_; n++) {
//      PsiZeta_[n] = PsiZeta_[n - 1]*(static_cast<FloatType>(n)*z_inv - D1[n -
//      1])
//                                   *(static_cast<FloatType>(n)*z_inv - D3[n -
//                                   1]);
//      D3[n] = D1[n] + std::complex<FloatType>(0.0, 1.0)/PsiZeta_[n];
//    }
//  }
//
//
//*******************************************************************************
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
//*******************************************************************************
//  template <typename FloatType>
//  void MultiLayerMie<FloatType>::calcPsiZeta(std::complex<FloatType> z,
//                                  std::vector<std::complex<FloatType> > &Psi,
//                                  std::vector<std::complex<FloatType> > &Zeta)
//                                  {
//
//    std::complex<FloatType> c_i(0.0, 1.0);
//    std::vector<std::complex<FloatType> > D1(nmax_ + 1), D3(nmax_ + 1);
//
//    // First, calculate the logarithmic derivatives
//    calcD1D3(z, D1, D3);
//
//    // Now, use the upward recurrence to calculate Psi and Zeta - equations
//    (20a) - (21b) Psi[0] = nmm::sin(z); Zeta[0] = nmm::sin(z) -
//    c_i*nmm::cos(z); for (int n = 1; n <= nmax_; n++) {
//      Psi[n]  =  Psi[n - 1]*(std::complex<FloatType>(n,0.0)/z - D1[n - 1]);
//      Zeta[n] = Zeta[n - 1]*(std::complex<FloatType>(n,0.0)/z - D3[n - 1]);
//    }
//  }

//**********************************************************************************
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
//**********************************************************************************
//  template <typename FloatType>
//  void MultiLayerMie<FloatType>::calcPiTau(const FloatType &costheta,
//                                std::vector<FloatType> &Pi,
//                                std::vector<FloatType> &Tau) {
//
//    int i;
//    //****************************************************//
//    // Equations (26a) - (26c)                            //
//    //****************************************************//
//    // Initialize Pi and Tau
//    Pi[0] = 1.0;  // n=1
//    Tau[0] = costheta;
//    // Calculate the actual values
//    if (nmax_ > 1) {
//      Pi[1] = 3*costheta*Pi[0]; //n=2
//      Tau[1] = 2*costheta*Pi[1] - 3*Pi[0];
//      for (i = 2; i < nmax_; i++) { //n=[3..nmax_]
//        Pi[i] = ((i + i + 1)*costheta*Pi[i - 1] - (i + 1)*Pi[i - 2])/i;
//        Tau[i] = (i + 1)*costheta*Pi[i] - (i + 2)*Pi[i - 1];
//      }
//    }
//  }  // end of MultiLayerMie::calcPiTau(...)

//**********************************************************************************
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
//**********************************************************************************
//  template <typename FloatType>
//  void MultiLayerMie<FloatType>::calcSpherHarm(const std::complex<FloatType>
//  Rho, const FloatType Theta, const FloatType Phi,
//                                    const std::complex<FloatType> &rn, const
//                                    std::complex<FloatType> &Dn, const
//                                    FloatType &Pi, const FloatType &Tau, const
//                                    FloatType &n,
//                                    std::vector<std::complex<FloatType> >&
//                                    Mo1n, std::vector<std::complex<FloatType>
//                                    >& Me1n,
//                                    std::vector<std::complex<FloatType> >&
//                                    No1n, std::vector<std::complex<FloatType>
//                                    >& Ne1n) {
//
//    // using eq 4.50 in BH
//    std::complex<FloatType> c_zero(0.0, 0.0);
//
//    using nmm::sin;
//    using nmm::cos;
//    Mo1n[0] = c_zero;
//    Mo1n[1] = cos(Phi)*Pi*rn/Rho;
//    Mo1n[2] = -sin(Phi)*Tau*rn/Rho;
//
//    Me1n[0] = c_zero;
//    Me1n[1] = -sin(Phi)*Pi*rn/Rho;
//    Me1n[2] = -cos(Phi)*Tau*rn/Rho;
//
//    No1n[0] = sin(Phi)*(n*n + n)*sin(Theta)*Pi*rn/Rho/Rho;
//    No1n[1] = sin(Phi)*Tau*Dn*rn/Rho;
//    No1n[2] = cos(Phi)*Pi*Dn*rn/Rho;
//
//    Ne1n[0] = cos(Phi)*(n*n + n)*sin(Theta)*Pi*rn/Rho/Rho;
//    Ne1n[1] = cos(Phi)*Tau*Dn*rn/Rho;
//    Ne1n[2] = -sin(Phi)*Pi*Dn*rn/Rho;
//  }  // end of MultiLayerMie::calcSpherHarm(...)
//

//******************************************************************************
// This functions calculate the logarithmic derivatives of the Riccati-Bessel
// functions (D1 and D3) for a complex argument (z).
// Equations (16a), (16b) and (18a) - (18d)
//
// Input parameters:
//   z: Complex argument to evaluate D1 and D3
//   nmax_: Maximum number of terms to calculate D1 and D3
//
// Helper to calculate nstar for scalar or SIMD
template <typename FloatType, typename Engine, typename ComplexType>
struct NStarCalculator {
  static int get(int nmax, ComplexType z, int valid_digits) {
    return nmie::getNStar<FloatType, Engine, ComplexType>(nmax, z, valid_digits);
  }
};

// Output parameters:
//   D1, D3: Logarithmic derivatives of the Riccati-Bessel functions
//******************************************************************************
template <typename FloatType, typename Engine = ScalarEngine<FloatType>, typename ComplexType = std::complex<FloatType>, typename ContainerType = std::vector<std::complex<FloatType> > >
void evalDownwardD1(const ComplexType z,
                    ContainerType& D1) {
  int lanes = Engine::Lanes();
  int nmax = (D1.size() / lanes) - 1;
  int valid_digits = 16;
#ifdef MULTI_PRECISION
  valid_digits += MULTI_PRECISION;
#endif
  int nstar = NStarCalculator<FloatType, Engine, ComplexType>::get(nmax, z, valid_digits);
  D1.resize((nstar + 1) * lanes);
  // Downward recurrence for D1 - equations (16a) and (16b)
  auto zero = Engine::make_complex(Engine::set(0.0), Engine::set(0.0));
  Engine::store(zero, &D1[nstar * lanes]);

  auto c_one = Engine::make_complex(Engine::set(1.0), Engine::set(0.0));
  auto z_inv = c_one / z;
  
  for (unsigned int n = nstar; n > 0; n--) {
    auto n_val = Engine::set(static_cast<FloatType>(n));
    auto n_complex = Engine::make_complex(n_val, Engine::set(0.0));
    auto term1 = n_complex * z_inv;
    
    auto d1_n = Engine::load(&D1[n * lanes]);
    auto denom = d1_n + term1;
    auto term2 = c_one / denom;
    
    auto res = term1 - term2;
    Engine::store(res, &D1[(n - 1) * lanes]);
  }
  // Use D1[0] from upward recurrence
  auto cot_z = complex_cot<FloatType, Engine>(z);
  Engine::store(cot_z, &D1[0]);
  D1.resize((nmax + 1) * lanes);
  //  printf("D1[0] = (%16.15g, %16.15g) z=(%16.15g,%16.15g)\n",
  //  D1[0].real(),D1[0].imag(),
  //         z.real(),z.imag());
}

//******************************************************************************
template <typename FloatType, typename Engine = ScalarEngine<FloatType>, typename ComplexType = std::complex<FloatType>, typename ContainerType = std::vector<std::complex<FloatType> > >
void evalUpwardD3(const ComplexType z,
                  const ContainerType& D1,
                  ContainerType& D3,
                  ContainerType& PsiZeta) {
  int lanes = Engine::Lanes();
  int nmax = (D1.size() / lanes) - 1;
  D3.resize((nmax + 1) * lanes);
  PsiZeta.resize((nmax + 1) * lanes);

  // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
  
  auto two = Engine::set(2.0);
  auto half = Engine::set(0.5);
  auto one = Engine::set(1.0);
  auto zero = Engine::set(0.0);
  
  auto z_re = Engine::get_real(z);
  auto z_im = Engine::get_imag(z);
  
  auto two_z_re = two * z_re;
  auto minus_two_z_im = Engine::set(-2.0) * z_im;
  
  auto cos_val = Engine::cos(two_z_re);
  auto sin_val = Engine::sin(two_z_re);
  auto exp_val = Engine::exp(minus_two_z_im);
  
  auto term_complex = Engine::make_complex(cos_val, sin_val);
  auto exp_complex = Engine::make_complex(exp_val, zero);
  
  auto prod = term_complex * exp_complex;
  
  auto one_complex = Engine::make_complex(one, zero);
  auto diff = one_complex - prod;
  
  auto half_complex = Engine::make_complex(half, zero);
  auto psi_zeta_0 = half_complex * diff;
  Engine::store(psi_zeta_0, &PsiZeta[0]);

  auto d3_0 = Engine::make_complex(zero, one);
  Engine::store(d3_0, &D3[0]);
  
  auto z_inv = one_complex / z;
  auto i_complex = Engine::make_complex(zero, one);
  
  for (int n = 1; n <= nmax; n++) {
    auto n_val = Engine::set(static_cast<FloatType>(n));
    auto n_complex = Engine::make_complex(n_val, zero);
    
    auto term = n_complex * z_inv;
    
    auto d1_nm1 = Engine::load(&D1[(n - 1) * lanes]);
    auto d3_nm1 = Engine::load(&D3[(n - 1) * lanes]);

    auto factor1 = term - d1_nm1;
    auto factor2 = term - d3_nm1;
    
    auto psi_zeta_nm1 = Engine::load(&PsiZeta[(n - 1) * lanes]);
    auto prod1 = psi_zeta_nm1 * factor1;
    auto psi_zeta_n = prod1 * factor2;
    Engine::store(psi_zeta_n, &PsiZeta[n * lanes]);
    
    auto term_div = i_complex / psi_zeta_n;
    auto d1_n = Engine::load(&D1[n * lanes]);
    auto d3_n = d1_n + term_div;
    Engine::store(d3_n, &D3[n * lanes]);
  }
}

template <typename FloatType, typename Engine = ScalarEngine<FloatType>, typename ComplexType = std::complex<FloatType> >
ComplexType complex_sin(const ComplexType z) {
  auto a = Engine::get_real(z);
  auto b = Engine::get_imag(z);
  
  auto zero = Engine::set(0.0);
  auto two = Engine::set(2.0);
  
  // term1 = (cos(a) + i*sin(a)) * exp(-b)
  auto cos_a = Engine::cos(a);
  auto sin_a = Engine::sin(a);
  auto exp_minus_b = Engine::exp(zero - b);
  
  auto e_ia = Engine::make_complex(cos_a, sin_a);
  auto term1 = e_ia * Engine::make_complex(exp_minus_b, zero);
  
  // term2 = (cos(-a) + i*sin(-a)) * exp(b)
  auto minus_a = zero - a;
  auto cos_minus_a = Engine::cos(minus_a);
  auto sin_minus_a = Engine::sin(minus_a);
  auto exp_b = Engine::exp(b);
  
  auto e_minus_ia = Engine::make_complex(cos_minus_a, sin_minus_a);
  auto term2 = e_minus_ia * Engine::make_complex(exp_b, zero);
  
  auto diff = term1 - term2;
  auto two_i = Engine::make_complex(zero, two);
  
  return diff / two_i;
}

//******************************************************************************
template <typename FloatType, typename Engine = ScalarEngine<FloatType>, typename ComplexType = std::complex<FloatType>, typename ContainerType = std::vector<std::complex<FloatType> > >
void evalUpwardPsi(const ComplexType z,
                   const ContainerType& D1,
                   ContainerType& Psi) {
  int lanes = Engine::Lanes();
  int nmax = (D1.size() / lanes) - 1;
  Psi.resize((nmax + 1) * lanes);

  // Now, use the upward recurrence to calculate Psi and Zeta - equations (20a)
  // - (21b)
  auto psi_0 = complex_sin<FloatType, Engine, ComplexType>(z);
  Engine::store(psi_0, &Psi[0]);
  
  auto one = Engine::set(1.0);
  auto zero = Engine::set(0.0);
  auto one_complex = Engine::make_complex(one, zero);
  auto z_inv = one_complex / z;
  
  for (int n = 1; n <= nmax; n++) {
    auto n_val = Engine::set(static_cast<FloatType>(n));
    auto n_complex = Engine::make_complex(n_val, zero);
    
    auto term = n_complex * z_inv;
    auto d1_nm1 = Engine::load(&D1[(n - 1) * lanes]);
    auto factor = term - d1_nm1;
    
    auto psi_nm1 = Engine::load(&Psi[(n - 1) * lanes]);
    auto psi_n = psi_nm1 * factor;
    Engine::store(psi_n, &Psi[n * lanes]);
  }
}

//******************************************************************************
template <typename FloatType, typename Engine = ScalarEngine<FloatType>, typename ComplexType = std::complex<FloatType> >
ComplexType complex_cos(const ComplexType z) {
  auto a = Engine::get_real(z);
  auto b = Engine::get_imag(z);
  
  auto zero = Engine::set(0.0);
  auto two = Engine::set(2.0);
  
  auto cos_a = Engine::cos(a);
  auto sin_a = Engine::sin(a);
  auto exp_minus_b = Engine::exp(zero - b);
  
  auto e_ia = Engine::make_complex(cos_a, sin_a);
  auto term1 = e_ia * Engine::make_complex(exp_minus_b, zero);
  
  auto minus_a = zero - a;
  auto cos_minus_a = Engine::cos(minus_a);
  auto sin_minus_a = Engine::sin(minus_a);
  auto exp_b = Engine::exp(b);
  
  auto e_minus_ia = Engine::make_complex(cos_minus_a, sin_minus_a);
  auto term2 = e_minus_ia * Engine::make_complex(exp_b, zero);
  
  auto sum = term1 + term2;
  auto two_complex = Engine::make_complex(two, zero);
  
  return sum / two_complex;
}

//******************************************************************************
// Sometimes in literature Zeta is also denoted as Ksi, it is a Riccati-Bessel
// function of third kind.
//******************************************************************************
template <typename FloatType, typename Engine = ScalarEngine<FloatType>, typename ComplexType = std::complex<FloatType>, typename ContainerType = std::vector<std::complex<FloatType> > >
void evalUpwardZeta(const ComplexType z,
                    const ContainerType& D3,
                    ContainerType& Zeta) {
  int nmax = Zeta.size() - 1;
  
  auto zero = Engine::set(0.0);
  auto one = Engine::set(1.0);
  auto i_complex = Engine::make_complex(zero, one);
  
  // Zeta[0] = sin(z) - i * cos(z)
  auto sin_z = complex_sin<FloatType, Engine, ComplexType>(z);
  auto cos_z = complex_cos<FloatType, Engine, ComplexType>(z);
  
  Zeta[0] = sin_z - i_complex * cos_z;
  
  auto one_complex = Engine::make_complex(one, zero);
  auto z_inv = one_complex / z;

  for (int n = 1; n <= nmax; n++) {
    auto n_val = Engine::set(static_cast<FloatType>(n));
    auto n_complex = Engine::make_complex(n_val, zero);
    
    auto factor = n_complex * z_inv - D3[n - 1];
    
    Zeta[n] = Zeta[n - 1] * factor;
  }
}

//******************************************************************************
// void evalForwardRiccatiBessel(const FloatType x, const FloatType first, const
// FloatType second,
//                               std::vector<FloatType> &values) {
//   values[0] = first;
//   values[1] = second;
//   int nmax = values.size();
//   for (int i = 1; i < nmax-1; i++) {
//     values[i+1] = (1 + 2*i) * values[i]/x - values[i-1];
//   }
// }
//
//******************************************************************************
// void evalChi(const FloatType x, std::vector<FloatType> &Chi) {
//   auto first = nmm::cos(x);
//   auto second = first/x + nmm::sin(x);
//   evalForwardRiccatiBessel(x, first, second, Chi);
// }
//
//******************************************************************************
// void evalPsi(const FloatType x, std::vector<FloatType> &Psi) {
//   auto first = nmm::sin(x);
//   auto second = first/x - nmm::cos(x);
//   evalForwardRiccatiBessel(x, first, second, Psi);
// }
//
//******************************************************************************
// void composeZeta(const std::vector<FloatType> &Psi,
//                  const std::vector<FloatType> &Chi,
//                  std::vector< std::complex<FloatType>> &Zeta) {
//   int nmax = Zeta.size();
//   for (int i = 0; i < nmax; i++) {
//     Zeta[i] = std::complex<FloatType > (Psi[i], Chi[i]);
//   }
// }
//******************************************************************************
template <typename FloatType, typename Engine = ScalarEngine<FloatType>>
void evalPsiZetaD1D3(const std::complex<FloatType> cxd,
                     std::vector<std::complex<FloatType>>& Psi,
                     std::vector<std::complex<FloatType>>& Zeta,
                     std::vector<std::complex<FloatType>>& D1,
                     std::vector<std::complex<FloatType>>& D3) {
  int nmax = Psi.size() - 1;
  std::vector<std::complex<FloatType>> PsiZeta(nmax + 1);

  for (int n = 0; n <= nmax; n++) {
    D1[n] = std::complex<FloatType>(0.0, -1.0);
    D3[n] = std::complex<FloatType>(0.0, 1.0);
  }

  evalDownwardD1<FloatType, Engine>(cxd, D1);
  evalUpwardPsi<FloatType, Engine>(cxd, D1, Psi);
  evalUpwardD3<FloatType, Engine>(cxd, D1, D3, PsiZeta);
  for (unsigned int i = 0; i < Zeta.size(); i++) {
    Zeta[i] = PsiZeta[i] / Psi[i];
  }
}

}  // end of namespace nmie
#endif  // SRC_SPECIAL_FUNCTIONS_IMPL_HPP_
