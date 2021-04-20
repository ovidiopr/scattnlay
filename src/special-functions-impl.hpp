#ifndef SRC_SPECIAL_FUNCTIONS_IMPL_HPP_
#define SRC_SPECIAL_FUNCTIONS_IMPL_HPP_
//**********************************************************************************//
//    Copyright (C) 2009-2018  Ovidio Pena <ovidio@bytesfall.com>                   //
//    Copyright (C) 2013-2018  Konstantin Ladutenko <kostyfisik@gmail.com>          //
//                                                                                  //
//    This file is part of scattnlay                                                //
//                                                                                  //
//    This program is free software: you can redistribute it and/or modify          //
//    it under the terms of the GNU General Public License as published by          //
//    the Free Software Foundation, either version 3 of the License, or             //
//    (at your option) any later version.                                           //
//                                                                                  //
//    This program is distributed in the hope that it will be useful,               //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of                //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                 //
//    GNU General Public License for more details.                                  //
//                                                                                  //
//    The only additional remark is that we expect that all publications            //
//    describing work using this software, or all commercial products               //
//    using it, cite at least one of the following references:                      //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
//    [2] K. Ladutenko, U. Pal, A. Rivera, and O. Pena-Rodriguez, "Mie              //
//        calculation of electromagnetic near-field for a multilayered              //
//        sphere," Computer Physics Communications, vol. 214, May 2017,             //
//        pp. 225-230.                                                              //
//                                                                                  //
//    You should have received a copy of the GNU General Public License             //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.         //
//**********************************************************************************//

//**********************************************************************************//
// This class implements the algorithm for a multilayered sphere described by:      //
//    [1] W. Yang, "Improved recursive algorithm for light scattering by a          //
//        multilayered sphere,” Applied Optics, vol. 42, Mar. 2003, pp. 1710-1720.  //
//                                                                                  //
// You can find the description of all the used equations in:                       //
//    [2] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
//    [3] K. Ladutenko, U. Pal, A. Rivera, and O. Pena-Rodriguez, "Mie              //
//        calculation of electromagnetic near-field for a multilayered              //
//        sphere," Computer Physics Communications, vol. 214, May 2017,             //
//        pp. 225-230.                                                              //
//                                                                                  //
// Hereinafter all equations numbers refer to [2]                                   //
//**********************************************************************************//
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>

#include "nmie.hpp"
#include "nmie-precision.hpp"

namespace nmie {
////helper functions
//template<class T> inline T pow2(const T value) {return value*value;}

// Note, that Kapteyn seems to be too optimistic (at least by 3 digits
// in some cases) for forward recurrence, see D1test with WYang_data
int evalKapteynNumberOfLostSignificantDigits(const int ni,
                                             const std::complex<FloatType> z) {
  using nmm::abs, nmm::imag, nmm::real, nmm::log, nmm::sqrt, nmm::round;
  auto n = static_cast<FloatType>(ni);
  auto one = std::complex<FloatType> (1, 0);
  return round((
      abs(imag(z)) - log(2) - n * ( real( log( z/n) + sqrt(one
      - pow2(z/n)) - log (one + sqrt (one
      - pow2(z/n)))
      )))/ log(10));
}

int getNStar(int nmax, std::complex<FloatType> z, const int valid_digits) {
  int nstar = nmax;
  int forwardLoss = evalKapteynNumberOfLostSignificantDigits(nmax, z);
  int increment = static_cast<int>(std::ceil(
      std::max(4* std::pow(std::abs(z), 1/3.0), 5.0)
      ));
  int backwardLoss =evalKapteynNumberOfLostSignificantDigits(nstar, z);
  while ( backwardLoss - forwardLoss < valid_digits) {
    nstar += increment;
    backwardLoss = evalKapteynNumberOfLostSignificantDigits(nstar,z);
  };
  return nstar;
}

// Custom implementation of complex cot function to avoid overflow
// if Im(z) < 0, then it evaluates cot(z) as conj(cot(conj(z)))
const std::complex<FloatType>
complex_cot(const std::complex<FloatType> z) {
  auto Remx = z.real();
  auto Immx = z.imag();
  int sign = (Immx>0) ? 1: -1; // use complex conj if needed for exp and return
  auto exp = nmm::exp(-2*sign*Immx);
  auto tan = nmm::tan(Remx);
  auto a = tan -  exp*tan;
  auto b =  1 + exp;
  auto c = -1 + exp;
  auto d = tan + exp* tan;
  auto c_one = std::complex<FloatType>(0, 1);
  return (a*c + b*d)/(pow2(c) + pow2(d)) +
      c_one * (sign* (b*c - a*d)/(pow2(c) + pow2(d)));
}

// Forward iteration for evaluation of ratio of the Riccati–Bessel functions
void evalForwardR (const std::complex<FloatType> z,
                   std::vector<std::complex<FloatType> >& r) {
  if (r.size() < 1) throw std::invalid_argument(
      "We need non-zero evaluations of ratio of the Riccati–Bessel functions.\n");
  // r0 = cot(z)
//  r[0] = nmm::cos(z)/nmm::sin(z);
  r[0] = complex_cot(z);
  for (unsigned int n = 0; n < r.size() - 1; n++) {
    r[n+1] = static_cast<FloatType>(1)/( static_cast<FloatType>(2*n+1)/z - r[n] );
  }
}


// Backward iteration for evaluation of ratio of the Riccati–Bessel functions
void evalBackwardR (const std::complex<FloatType> z,
                   std::vector<std::complex<FloatType> >& r) {
  if (r.size() < 1) throw std::invalid_argument(
        "We need non-zero evaluations of ratio of the Riccati–Bessel functions.\n");
  int nmax = r.size()-1 ;
  auto tmp = static_cast<FloatType>(2*nmax + 1)/z;
  r[nmax] = tmp;
  for (int n = nmax -1; n >= 0; n--) {
    r[n] = -static_cast<FloatType>(1)/r[n+1] + static_cast<FloatType>(2*n+1)/z;
  }
//  r[0] = nmm::cos(z)/nmm::sin(z);
//  nmm::cout << "R0 = " << r[0] <<" at arg = "<<z<<'\n';
}

void convertRtoD1(const std::complex<FloatType> z,
                  std::vector<std::complex<FloatType> >& r,
                  std::vector<std::complex<FloatType> >& D1) {
  if (D1.size() > r.size()) throw std::invalid_argument(
      "Not enough elements in array of ratio of the Riccati–Bessel functions to convert it into logarithmic derivative array.\n");
  std::complex<FloatType> Dold;
  for (unsigned int n = 0; n < D1.size(); n++) {
    Dold = D1[n];
    D1[n] = r[n] - static_cast<FloatType>(n)/z;
//    std::cout << "D1old - D1[n] :" << Dold-D1[n] << '\n';
  }
}

// ********************************************************************** //
void evalForwardD (const std::complex<FloatType> z,
                     std::vector<std::complex<FloatType> >& D) {
  int nmax = D.size();
  FloatType one = static_cast<FloatType >(1);
  for (int n = 1; n < nmax; n++) {
    FloatType nf = static_cast<FloatType > (n);
    D[n] = one/ (nf/z  - D[n-1]) - nf/z;
  }
}

// ********************************************************************** //
void evalForwardD1 (const std::complex<FloatType> z,
                   std::vector<std::complex<FloatType> >& D) {
  if (D.size()<1) throw std::invalid_argument("Should have a leas one element!\n");
  D[0] = std::cos(z)/std::sin(z);
  evalForwardD(z,D);
}

//  template <typename FloatType>
//  void MultiLayerMie<FloatType>::calcD1D3(const std::complex<FloatType> z,
//                               std::vector<std::complex<FloatType> >& D1,
//                               std::vector<std::complex<FloatType> >& D3) {
//
//    // if (cabs(D1[0]) > 1.0e15) {
//    //   throw std::invalid_argument("Unstable D1! Please, try to change input parameters!\n");
//    // //printf("Warning: Potentially unstable D1! Please, try to change input parameters!\n");
//    // }
//
//    // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
//    PsiZeta_[0] = static_cast<FloatType>(0.5)*(static_cast<FloatType>(1.0) - std::complex<FloatType>(nmm::cos(2.0*z.real()), nmm::sin(2.0*z.real()))
//                 *static_cast<FloatType>(nmm::exp(-2.0*z.imag())));
//    D3[0] = std::complex<FloatType>(0.0, 1.0);
//
//    for (int n = 1; n <= nmax_; n++) {
//      PsiZeta_[n] = PsiZeta_[n - 1]*(static_cast<FloatType>(n)*zinv - D1[n - 1])
//                                   *(static_cast<FloatType>(n)*zinv - D3[n - 1]);
//      D3[n] = D1[n] + std::complex<FloatType>(0.0, 1.0)/PsiZeta_[n];
//    }
//  }
//
//
  //**********************************************************************************//
  // This function calculates the Riccati-Bessel functions (Psi and Zeta) for a       //
  // complex argument (z).                                                            //
  // Equations (20a) - (21b)                                                          //
  //                                                                                  //
  // Input parameters:                                                                //
  //   z: Complex argument to evaluate Psi and Zeta                                   //
  //   nmax: Maximum number of terms to calculate Psi and Zeta                        //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Psi, Zeta: Riccati-Bessel functions                                            //
  //**********************************************************************************//
//  template <typename FloatType>
//  void MultiLayerMie<FloatType>::calcPsiZeta(std::complex<FloatType> z,
//                                  std::vector<std::complex<FloatType> >& Psi,
//                                  std::vector<std::complex<FloatType> >& Zeta) {
//
//    std::complex<FloatType> c_i(0.0, 1.0);
//    std::vector<std::complex<FloatType> > D1(nmax_ + 1), D3(nmax_ + 1);
//
//    // First, calculate the logarithmic derivatives
//    calcD1D3(z, D1, D3);
//
//    // Now, use the upward recurrence to calculate Psi and Zeta - equations (20a) - (21b)
//    Psi[0] = std::sin(z);
//    Zeta[0] = std::sin(z) - c_i*std::cos(z);
//    for (int n = 1; n <= nmax_; n++) {
//      Psi[n]  =  Psi[n - 1]*(std::complex<FloatType>(n,0.0)/z - D1[n - 1]);
//      Zeta[n] = Zeta[n - 1]*(std::complex<FloatType>(n,0.0)/z - D3[n - 1]);
//    }
//  }


  //**********************************************************************************//
  // This function calculates Pi and Tau for a given value of cos(Theta).             //
  // Equations (26a) - (26c)                                                          //
  //                                                                                  //
  // Input parameters:                                                                //
  //   nmax_: Maximum number of terms to calculate Pi and Tau                         //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Pi, Tau: Angular functions Pi and Tau, as defined in equations (26a) - (26c)   //
  //**********************************************************************************//
//  template <typename FloatType>
//  void MultiLayerMie<FloatType>::calcPiTau(const FloatType& costheta,
//                                std::vector<FloatType>& Pi, std::vector<FloatType>& Tau) {
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


  //**********************************************************************************//
  // This function calculates vector spherical harmonics (eq. 4.50, p. 95 BH),        //
  // required to calculate the near-field parameters.                                 //
  //                                                                                  //
  // Input parameters:                                                                //
  //   Rho: Radial distance                                                           //
  //   Phi: Azimuthal angle                                                           //
  //   Theta: Polar angle                                                             //
  //   rn: Either the spherical Ricatti-Bessel function of first or third kind        //
  //   Dn: Logarithmic derivative of rn                                               //
  //   Pi, Tau: Angular functions Pi and Tau                                          //
  //   n: Order of vector spherical harmonics                                         //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Mo1n, Me1n, No1n, Ne1n: Complex vector spherical harmonics                     //
  //**********************************************************************************//
//  template <typename FloatType>
//  void MultiLayerMie<FloatType>::calcSpherHarm(const std::complex<FloatType> Rho, const FloatType Theta, const FloatType Phi,
//                                    const std::complex<FloatType>& rn, const std::complex<FloatType>& Dn,
//                                    const FloatType& Pi, const FloatType& Tau, const FloatType& n,
//                                    std::vector<std::complex<FloatType> >& Mo1n, std::vector<std::complex<FloatType> >& Me1n,
//                                    std::vector<std::complex<FloatType> >& No1n, std::vector<std::complex<FloatType> >& Ne1n) {
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

//**********************************************************************************//
// This functions calculate the logarithmic derivatives of the Riccati-Bessel       //
// functions (D1 and D3) for a complex argument (z).                                //
// Equations (16a), (16b) and (18a) - (18d)                                         //
//                                                                                  //
// Input parameters:                                                                //
//   z: Complex argument to evaluate D1 and D3                                      //
//   nmax_: Maximum number of terms to calculate D1 and D3                          //
//                                                                                  //
// Output parameters:                                                               //
//   D1, D3: Logarithmic derivatives of the Riccati-Bessel functions                //
//**********************************************************************************//
void evalDownwardD1 (const std::complex<FloatType> z,
                     std::vector<std::complex<FloatType> >& D1) {
  int nmax = D1.size() - 1;
  // Downward recurrence for D1 - equations (16a) and (16b)
  D1[nmax] = std::complex<FloatType>(0.0, 0.0);
  std::complex<FloatType> c_one(1.0, 0.0);
  const std::complex<FloatType> zinv = std::complex<FloatType>(1.0, 0.0)/z;
  for (unsigned int n = nmax; n > 0; n--) {
    D1[n - 1] = static_cast<FloatType>(n)*zinv - c_one/
        (D1[n] + static_cast<FloatType>(n)*zinv);
  }
  // Use D1[0] from upward recurrence
  D1[0] = complex_cot(z);

//  printf("D1[0] = (%16.15g, %16.15g) z=(%16.15g,%16.15g)\n", D1[0].real(),D1[0].imag(),
//         z.real(),z.imag());
}


void evalUpwardD3 (const std::complex<FloatType> z,
                   const std::vector<std::complex<FloatType> >& D1,
                   std::vector<std::complex<FloatType> >& D3,
                   std::vector<std::complex<FloatType> >& PsiZeta) {
  int nmax = D1.size()-1;
  // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
  PsiZeta[0] = static_cast<FloatType>(0.5)*(static_cast<FloatType>(1.0) - std::complex<FloatType>(nmm::cos(2.0*z.real()), nmm::sin(2.0*z.real()))
      *static_cast<FloatType>(nmm::exp(-2.0*z.imag())));
  D3[0] = std::complex<FloatType>(0.0, 1.0);
  const std::complex<FloatType> zinv = std::complex<FloatType>(1.0, 0.0)/z;
  for (int n = 1; n <= nmax; n++) {
    PsiZeta[n] = PsiZeta[n - 1]*(static_cast<FloatType>(n)*zinv - D1[n - 1])
        *(static_cast<FloatType>(n)*zinv - D3[n - 1]);
    D3[n] = D1[n] + std::complex<FloatType>(0.0, 1.0)/PsiZeta[n];
  }
}


void evalUpwardPsi (const std::complex<FloatType> z,
                    const std::vector<std::complex<FloatType> > D1,
                   std::vector<std::complex<FloatType> >& Psi) {
  int nmax = Psi.size() - 1;
  std::complex<FloatType> c_i(0.0, 1.0);
  // Now, use the upward recurrence to calculate Psi and Zeta - equations (20a) - (21b)
  Psi[0] = std::sin(z);
  for (int n = 1; n <= nmax; n++) {
    Psi[n]  =  Psi[n - 1]*(std::complex<FloatType>(n,0.0)/z - D1[n - 1]);
  }
}

// Sometimes in literature Zeta is also denoted as Ksi, it is a Riccati-Bessel function of third kind.
void evalUpwardZeta (const std::complex<FloatType> z,
                    const std::vector<std::complex<FloatType> > D3,
                    std::vector<std::complex<FloatType> >& Zeta) {
  int nmax = Zeta.size() - 1;
  std::complex<FloatType> c_i(0.0, 1.0);
  // Now, use the upward recurrence to calculate Zeta and Zeta - equations (20a) - (21b)
  Zeta[0] = std::sin(z) - c_i*std::cos(z);
  for (int n = 1; n <= nmax; n++) {
    Zeta[n]  =  Zeta[n - 1]*(std::complex<FloatType>(n, 0.0)/z - D3[n - 1]);
  }
}

//void evalForwardRiccatiBessel(const FloatType x, const FloatType first, const FloatType second,
//                              std::vector<FloatType> &values) {
//  values[0] = first;
//  values[1] = second;
//  int nmax = values.size();
//  for (int i = 1; i < nmax-1; i++) {
//    values[i+1] = (1 + 2*i) * values[i]/x - values[i-1];
//  }
//}
//
//void evalChi(const FloatType x, std::vector<FloatType> &Chi) {
//  auto first = nmm::cos(x);
//  auto second = first/x + nmm::sin(x);
//  evalForwardRiccatiBessel(x, first, second, Chi);
//}
//
//void evalPsi(const FloatType x, std::vector<FloatType> &Psi) {
//  auto first = nmm::sin(x);
//  auto second = first/x - nmm::cos(x);
//  evalForwardRiccatiBessel(x, first, second, Psi);
//}
//
//void composeZeta(const std::vector<FloatType> &Psi,
//                 const std::vector<FloatType> &Chi,
//                 std::vector< std::complex<FloatType>> &Zeta) {
//  int nmax = Zeta.size();
//  for (int i = 0; i < nmax; i++) {
//    Zeta[i] = std::complex<FloatType > (Psi[i], Chi[i]);
//  }
//}

}  // end of namespace nmie
#endif  // SRC_SPECIAL_FUNCTIONS_IMPL_HPP_
