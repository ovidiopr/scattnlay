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

template <typename FloatType, typename Engine = ScalarEngine, typename ComplexType>
ComplexType calc_an(int n,
                      typename Engine::RealV XL, 
                      ComplexType Ha,
                      ComplexType mL,
                      ComplexType PsiXL,
                      ComplexType ZetaXL,
                      ComplexType PsiXLM1,
                      ComplexType ZetaXLM1) {
    auto n_val = Engine::set(static_cast<FloatType>(n));
    auto zero = Engine::set(0.0);
    auto n_complex = Engine::make_complex(n_val, zero);
    auto XL_complex = Engine::make_complex(XL, zero);
    
    // (Ha / mL + n / XL)
    auto term1 = Engine::add(Engine::div(Ha, mL), Engine::div(n_complex, XL_complex));
    
    // Num = term1 * PsiXL - PsiXLM1
    auto Num = Engine::sub(Engine::mul(term1, PsiXL), PsiXLM1);
    
    // Denom = term1 * ZetaXL - ZetaXLM1
    auto Denom = Engine::sub(Engine::mul(term1, ZetaXL), ZetaXLM1);

    return Engine::div(Num, Denom);
}

template <typename FloatType, typename Engine = ScalarEngine, typename ComplexType>
ComplexType calc_bn(int n,
                      typename Engine::RealV XL, 
                      ComplexType Hb,
                      ComplexType mL,
                      ComplexType PsiXL,
                      ComplexType ZetaXL,
                      ComplexType PsiXLM1,
                      ComplexType ZetaXLM1) {
    auto n_val = Engine::set(static_cast<FloatType>(n));
    auto zero = Engine::set(0.0);
    auto n_complex = Engine::make_complex(n_val, zero);
    auto XL_complex = Engine::make_complex(XL, zero);
    
    // (mL * Hb + n / XL)
    auto term1 = Engine::add(Engine::mul(mL, Hb), Engine::div(n_complex, XL_complex));
    
    // Num = term1 * PsiXL - PsiXLM1
    auto Num = Engine::sub(Engine::mul(term1, PsiXL), PsiXLM1);
    
    // Denom = term1 * ZetaXL - ZetaXLM1
    auto Denom = Engine::sub(Engine::mul(term1, ZetaXL), ZetaXLM1);

    return Engine::div(Num, Denom);
}

// class implementation

// ********************************************************************** //
// Returns previously calculated Qext                                     //
// ********************************************************************** //
template <typename FloatType>
template <typename outputType>
outputType MultiLayerMie<FloatType>::GetQext() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return static_cast<outputType>(Qext_);
}

// ********************************************************************** //
// Returns previously calculated Qabs                                     //
// ********************************************************************** //
template <typename FloatType>
template <typename outputType>
outputType MultiLayerMie<FloatType>::GetQabs() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return static_cast<outputType>(Qabs_);
}

// ********************************************************************** //
// Returns previously calculated Qsca                                     //
// ********************************************************************** //
template <typename FloatType>
template <typename outputType>
outputType MultiLayerMie<FloatType>::GetQsca() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return static_cast<outputType>(Qsca_);
}

// ********************************************************************** //
// Returns previously calculated Qbk                                      //
// ********************************************************************** //
template <typename FloatType>
template <typename outputType>
outputType MultiLayerMie<FloatType>::GetQbk() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return static_cast<outputType>(Qbk_);
}

// ********************************************************************** //
// Returns previously calculated Qpr                                      //
// ********************************************************************** //
template <typename FloatType>
template <typename outputType>
outputType MultiLayerMie<FloatType>::GetQpr() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return static_cast<outputType>(Qpr_);
}

// ********************************************************************** //
// Returns previously calculated asymmetry factor                         //
// ********************************************************************** //
template <typename FloatType>
template <typename outputType>
outputType MultiLayerMie<FloatType>::GetAsymmetryFactor() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return static_cast<outputType>(asymmetry_factor_);
}

// ********************************************************************** //
// Returns previously calculated Albedo                                   //
// ********************************************************************** //
template <typename FloatType>
template <typename outputType>
outputType MultiLayerMie<FloatType>::GetAlbedo() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return static_cast<outputType>(albedo_);
}

// ********************************************************************** //
// Returns previously calculated S1                                       //
// ********************************************************************** //
template <typename FloatType>
std::vector<std::complex<FloatType>> MultiLayerMie<FloatType>::GetS1() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return S1_;
}

// ********************************************************************** //
// Returns previously calculated S2                                       //
// ********************************************************************** //
template <typename FloatType>
std::vector<std::complex<FloatType>> MultiLayerMie<FloatType>::GetS2() {
  if (!isMieCalculated_)
    throw std::invalid_argument(
        "You should run calculations before result request!");
  return S2_;
}

// ********************************************************************** //
// Modify scattering (theta) angles                                       //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMie<FloatType>::SetAngles(const std::vector<FloatType>& angles) {
  MarkUncalculated();
  theta_ = angles;
}

// ********************************************************************** //
// Modify size of all layers                                             //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMie<FloatType>::SetLayersSize(
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
template <typename FloatType>
void MultiLayerMie<FloatType>::SetLayersIndex(
    const std::vector<std::complex<FloatType>>& index) {
  MarkUncalculated();
  refractive_index_ = index;
}

// ********************************************************************** //
// Modify coordinates for field calculation                               //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMie<FloatType>::SetFieldCoords(
    const std::vector<std::vector<FloatType>>& coords) {
  if (coords.size() != 3)
    throw std::invalid_argument(
        "Error! Wrong dimension of field monitor points!");
  if (coords[0].size() != coords[1].size() ||
      coords[0].size() != coords[2].size())
    throw std::invalid_argument(
        "Error! Missing coordinates for field monitor points!");
  coords_ = coords;
}

// ********************************************************************** //
// Modify index of PEC layer                                              //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMie<FloatType>::SetPECLayer(int layer_position) {
  MarkUncalculated();
  if (layer_position < 0 && layer_position != -1)
    throw std::invalid_argument("Error! Layers are numbered from 0!");
  PEC_layer_position_ = layer_position;
}

// ********************************************************************** //
// Set maximun number of terms to be used                                 //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMie<FloatType>::SetMaxTerms(int nmax) {
  MarkUncalculated();
  nmax_preset_ = nmax;
}

// ********************************************************************** //
// Get total size parameter of particle                                   //
// ********************************************************************** //
template <typename FloatType>
FloatType MultiLayerMie<FloatType>::GetSizeParameter() {
  if (size_param_.size() > 0)
    return size_param_.back();
  else
    return 0;
}

// ********************************************************************** //
// Mark uncalculated                                                      //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMie<FloatType>::MarkUncalculated() {
  isExpCoeffsCalc_ = false;
  isScaCoeffsCalc_ = false;

  isMieCalculated_ = false;
}
// ********************************************************************** //
// Clear layer information                                                //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMie<FloatType>::ClearLayers() {
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
template <typename FloatType>
unsigned int MultiLayerMie<FloatType>::calcNstop(FloatType xL) {
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
template <typename FloatType>
unsigned int MultiLayerMie<FloatType>::calcNmax(FloatType xL) {
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
template <typename FloatType>
std::complex<FloatType> MultiLayerMie<FloatType>::calc_an(
    int n,
    FloatType XL,
    std::complex<FloatType> Ha,
    std::complex<FloatType> mL,
    std::complex<FloatType> PsiXL,
    std::complex<FloatType> ZetaXL,
    std::complex<FloatType> PsiXLM1,
    std::complex<FloatType> ZetaXLM1) {
  return nmie::calc_an<FloatType, ScalarEngine, std::complex<FloatType>>(
      n, XL, Ha, mL, PsiXL, ZetaXL, PsiXLM1, ZetaXLM1);
}

// ********************************************************************** //
// Calculate bn - equation (6)                                            //
// ********************************************************************** //
template <typename FloatType>
std::complex<FloatType> MultiLayerMie<FloatType>::calc_bn(
    int n,
    FloatType XL,
    std::complex<FloatType> Hb,
    std::complex<FloatType> mL,
    std::complex<FloatType> PsiXL,
    std::complex<FloatType> ZetaXL,
    std::complex<FloatType> PsiXLM1,
    std::complex<FloatType> ZetaXLM1) {
  return nmie::calc_bn<FloatType, ScalarEngine, std::complex<FloatType>>(
      n, XL, Hb, mL, PsiXL, ZetaXL, PsiXLM1, ZetaXLM1);
}

// ********************************************************************** //
// Calculates S1 - equation (25a)                                         //
// ********************************************************************** //
template <typename FloatType>
std::complex<FloatType> MultiLayerMie<FloatType>::calc_S1(
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
template <typename FloatType>
std::complex<FloatType> MultiLayerMie<FloatType>::calc_S2(
    int n,
    std::complex<FloatType> an,
    std::complex<FloatType> bn,
    FloatType Pi,
    FloatType Tau) {
  return calc_S1(n, an, bn, Tau, Pi);
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
void MultiLayerMie<FloatType>::calcD1D3(
    const std::complex<FloatType> z,
    std::vector<std::complex<FloatType>>& D1,
    std::vector<std::complex<FloatType>>& D3) {
  std::vector<std::complex<FloatType>> PsiZeta(nmax_ + 1);
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
void MultiLayerMie<FloatType>::calcPsiZeta(
    std::complex<FloatType> z,
    std::vector<std::complex<FloatType>>& Psi,
    std::vector<std::complex<FloatType>>& Zeta) {
  std::vector<std::complex<FloatType>> D1(nmax_ + 1), D3(nmax_ + 1),
      PsiZeta(nmax_ + 1);
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

template <typename FloatType>
void MultiLayerMie<FloatType>::calcPiTauAllTheta(
    const double from_Theta,
    const double to_Theta,
    std::vector<std::vector<FloatType>>& Pi,
    std::vector<std::vector<FloatType>>& Tau) {
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
template <typename FloatType>
void MultiLayerMie<FloatType>::calcPiTau(const FloatType& costheta,
                                         std::vector<FloatType>& Pi,
                                         std::vector<FloatType>& Tau) {
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
template <typename FloatType>
template <typename evalType>
void MultiLayerMie<FloatType>::calcSpherHarm(
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
    std::vector<std::complex<evalType>>& Ne1n) {
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
template <typename FloatType>
void MultiLayerMie<FloatType>::calcScattCoeffs() {
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
  int fl = (pl > 0) ? pl : 0;
  if (nmax_preset_ <= 0)
    nmax_ = calcNmax();
  else
    nmax_ = nmax_preset_;

  std::complex<FloatType> z1, z2;
  //**************************************************************************//
  // Note that since Fri, Nov 14, 2014 all arrays start from 0 (zero), which  //
  // means that index = layer number - 1 or index = n - 1. The only exception //
  // are the arrays for representing D1, D3 and Q because they need a value   //
  // for the index 0 (zero), hence it is important to consider this shift     //
  // between different arrays. The change was done to optimize memory usage.  //
  //**************************************************************************//
  // Allocate memory to the arrays
  std::vector<std::complex<FloatType>> D1_mlxl(nmax_ + 1), D1_mlxlM1(nmax_ + 1),
      D3_mlxl(nmax_ + 1), D3_mlxlM1(nmax_ + 1);

  std::vector<std::vector<std::complex<FloatType>>> Q(L), Ha(L), Hb(L);

  for (int l = 0; l < L; l++) {
    Q[l].resize(nmax_ + 1, static_cast<FloatType>(0.0));
    Ha[l].resize(nmax_, static_cast<FloatType>(0.0));
    Hb[l].resize(nmax_, static_cast<FloatType>(0.0));
  }

  an_.resize(nmax_, static_cast<FloatType>(0.0));
  bn_.resize(nmax_, static_cast<FloatType>(0.0));

  std::vector<std::complex<FloatType>> PsiXL(nmax_ + 1), ZetaXL(nmax_ + 1);

  //*************************************************//
  // Calculate D1 and D3 for z1 in the first layer   //
  //*************************************************//
  if (fl == pl) {  // PEC layer
    for (int n = 0; n <= nmax_; n++) {
      D1_mlxl[n] = std::complex<FloatType>(0.0, -1.0);
      D3_mlxl[n] = std::complex<FloatType>(0.0, 1.0);
    }
  } else {  // Regular layer
    z1 = x[fl] * m[fl];
    // Calculate D1 and D3
    calcD1D3(z1, D1_mlxl, D3_mlxl);
  }

  //******************************************************************//
  // Calculate Ha and Hb in the first layer - equations (7a) and (8a) //
  //******************************************************************//
  for (int n = 0; n < nmax_; n++) {
    Ha[fl][n] = D1_mlxl[n + 1];
    Hb[fl][n] = D1_mlxl[n + 1];
  }
  //*****************************************************//
  // Iteration from the second layer to the last one (L) //
  //*****************************************************//
  std::complex<FloatType> Temp, Num, Denom;
  std::complex<FloatType> G1, G2;
  for (int l = fl + 1; l < L; l++) {
    //************************************************************//
    // Calculate D1 and D3 for z1 and z2 in the layers fl + 1..L   //
    //************************************************************//
    z1 = x[l] * m[l];
    z2 = x[l - 1] * m[l];
    // Calculate D1 and D3 for z1
    calcD1D3(z1, D1_mlxl, D3_mlxl);
    // Calculate D1 and D3 for z2
    calcD1D3(z2, D1_mlxlM1, D3_mlxlM1);

    //*************************************************//
    // Calculate Q, Ha and Hb in the layers fl + 1..L   //
    //*************************************************//
    // Upward recurrence for Q - equations (19a) and (19b)
    Num =
        std::complex<FloatType>(nmm::exp(-2.0 * (z1.imag() - z2.imag())), 0.0) *
        std::complex<FloatType>(
            nmm::cos(-2.0 * z2.real()) - nmm::exp(-2.0 * z2.imag()),
            nmm::sin(-2.0 * z2.real()));
    Denom = std::complex<FloatType>(
        nmm::cos(-2.0 * z1.real()) - nmm::exp(-2.0 * z1.imag()),
        nmm::sin(-2.0 * z1.real()));
    Q[l][0] = Num / Denom;

    for (int n = 1; n <= nmax_; n++) {
      Num = (z1 * D1_mlxl[n] + FloatType(n)) *
            (FloatType(n) - z1 * D3_mlxl[n - 1]);
      Denom = (z2 * D1_mlxlM1[n] + FloatType(n)) *
              (FloatType(n) - z2 * D3_mlxlM1[n - 1]);
      Q[l][n] = ((pow2(x[l - 1] / x[l]) * Q[l][n - 1]) * Num) / Denom;
    }
    // Upward recurrence for Ha and Hb - equations (7b), (8b) and (12) - (15)
    for (int n = 1; n <= nmax_; n++) {
      // Ha
      if ((l - 1) == pl) {  // The layer below the current one is a PEC layer
        G1 = -D1_mlxlM1[n];
        G2 = -D3_mlxlM1[n];
      } else {
        G1 = (m[l] * Ha[l - 1][n - 1]) - (m[l - 1] * D1_mlxlM1[n]);
        G2 = (m[l] * Ha[l - 1][n - 1]) - (m[l - 1] * D3_mlxlM1[n]);
      }  // end of if PEC
      Temp = Q[l][n] * G1;
      Num = (G2 * D1_mlxl[n]) - (Temp * D3_mlxl[n]);
      Denom = G2 - Temp;
      Ha[l][n - 1] = Num / Denom;
      // Hb
      if ((l - 1) == pl) {  // The layer below the current one is a PEC layer
        G1 = Hb[l - 1][n - 1];
        G2 = Hb[l - 1][n - 1];
      } else {
        G1 = (m[l - 1] * Hb[l - 1][n - 1]) - (m[l] * D1_mlxlM1[n]);
        G2 = (m[l - 1] * Hb[l - 1][n - 1]) - (m[l] * D3_mlxlM1[n]);
      }  // end of if PEC

      Temp = Q[l][n] * G1;
      Num = (G2 * D1_mlxl[n]) - (Temp * D3_mlxl[n]);
      Denom = (G2 - Temp);
      Hb[l][n - 1] = (Num / Denom);
    }  // end of for Ha and Hb terms
  }    // end of for layers iteration

  //**************************************//
  // Calculate Psi and Zeta for XL         //
  //**************************************//
  // Calculate PsiXL and ZetaXL
  calcPsiZeta(std::complex<FloatType>(x[L - 1], 0.0), PsiXL, ZetaXL);

  //*********************************************************************//
  // Finally, we calculate the scattering coefficients (an and bn) and   //
  // the angular functions (Pi and Tau). Note that for these arrays the  //
  // first layer is 0 (zero), in future versions all arrays will follow  //
  // this convention to save memory. (13 Nov, 2014)                      //
  //*********************************************************************//
  FloatType a0 = 0, b0 = 0;
  for (int n = 0; n < nmax_; n++) {
    //********************************************************************//
    // Expressions for calculating an and bn coefficients are not valid if //
    // there is only one PEC layer (ie, for a simple PEC sphere).          //
    //********************************************************************//
    if (pl < (L - 1)) {
      an_[n] = calc_an(n + 1, x[L - 1], Ha[L - 1][n], m[L - 1], PsiXL[n + 1],
                       ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
      bn_[n] = calc_bn(n + 1, x[L - 1], Hb[L - 1][n], m[L - 1], PsiXL[n + 1],
                       ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
    } else {
      an_[n] = calc_an(n + 1, x[L - 1], std::complex<FloatType>(0.0, 0.0),
                       std::complex<FloatType>(1.0, 0.0), PsiXL[n + 1],
                       ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
      bn_[n] = PsiXL[n + 1] / ZetaXL[n + 1];
    }
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

    // TODO seems to provide not enough terms for near-field calclulation.
    //      if (cabs(an_[n]) / a0 < convergence_threshold_ &&
    //          cabs(bn_[n]) / b0 < convergence_threshold_) {
    //        if (nmax_preset_ <= 0) nmax_ = n;
    //        break;
    //      }

    if (nmm::isnan(an_[n].real()) || nmm::isnan(an_[n].imag()) ||
        nmm::isnan(bn_[n].real()) || nmm::isnan(bn_[n].imag())) {
      std::cout
          << "nmax value was changed due to unexpected error!!! New values is "
          << n << " (was " << nmax_ << ")" << std::endl;
      nmax_ = n;
      break;
    }

  }  // end of for an and bn terms
  isScaCoeffsCalc_ = true;
}  // end of MultiLayerMie::calcScattCoeffs()

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
template <typename FloatType>
void MultiLayerMie<FloatType>::RunMieCalculation() {
  if (size_param_.size() != refractive_index_.size())
    throw std::invalid_argument(
        "Each size parameter should have only one index!");
  if (size_param_.size() == 0)
    throw std::invalid_argument("Initialize model first!");

  const std::vector<FloatType>& x = size_param_;

  // MarkUncalculated();

  // Calculate scattering coefficients
  if (!isScaCoeffsCalc_)
    calcScattCoeffs();

  // Initialize the scattering parameters
  Qext_ = 0.0;
  Qsca_ = 0.0;
  Qabs_ = 0.0;
  Qbk_ = 0.0;
  Qpr_ = 0.0;

  asymmetry_factor_ = 0.0;
  albedo_ = 0.0;

  // Initialize the scattering amplitudes
  std::vector<std::complex<FloatType>> tmp1(theta_.size(),
                                            std::complex<FloatType>(0.0, 0.0));
  S1_.swap(tmp1);
  S2_ = S1_;
  // Precalculate cos(theta) - gives about 5% speed up.
  std::vector<FloatType> costheta(theta_.size(), static_cast<FloatType>(0.0));
  for (unsigned int t = 0; t < theta_.size(); t++) {
    costheta[t] = nmm::cos(theta_[t]);
  }

  std::vector<FloatType> Pi(nmax_), Tau(nmax_);

  std::complex<FloatType> Qbktmp(0.0, 0.0);
  std::vector<std::complex<FloatType>> Qbktmp_ch(nmax_ - 1, Qbktmp);
  // By using downward recurrence we avoid loss of precision due to float
  // rounding errors See:
  // https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
  //      http://en.wikipedia.org/wiki/Loss_of_significance
  for (int n = nmax_ - 2; n >= 0; n--) {
    //      for (int n = 0; n < nmax_; n++) {
    const int n1 = n + 1;
    if (mode_n_ == Modes::kAll) {
      // Equation (27)
      Qext_ += (n1 + n1 + 1.0) * (an_[n].real() + bn_[n].real());
      // Equation (28)
      Qsca_ += (n1 + n1 + 1.0) *
               (an_[n].real() * an_[n].real() + an_[n].imag() * an_[n].imag() +
                bn_[n].real() * bn_[n].real() + bn_[n].imag() * bn_[n].imag());
      //        std::cout<<"n ="<< n1 << " ext:"<<Qext_ <<"
      //        sca:"<<Qsca_<<std::endl;
      // Equation (29)
      Qpr_ += ((n1 * (n1 + 2.0) / (n1 + 1.0)) *
                   ((an_[n] * std::conj(an_[n1]) + bn_[n] * std::conj(bn_[n1]))
                        .real()) +
               ((n1 + n1 + 1.0) / (n1 * (n1 + 1.0))) *
                   (an_[n] * std::conj(bn_[n])).real());
      // Equation (33)
      Qbktmp += (FloatType)(n1 + n1 + 1.0) * (1.0 - 2.0 * (n1 % 2)) *
                (an_[n] - bn_[n]);
      // Calculate the scattering amplitudes (S1 and S2) Equations (25a) - (25b)
      for (unsigned int t = 0; t < theta_.size(); t++) {
        calcPiTau(costheta[t], Pi, Tau);
        S1_[t] += calc_S1(n1, an_[n], bn_[n], Pi[n], Tau[n]);
        S2_[t] += calc_S2(n1, an_[n], bn_[n], Pi[n], Tau[n]);
      }
      continue;
    }
    if (n1 == mode_n_) {
      if (mode_type_ == Modes::kElectric || mode_type_ == Modes::kAll) {
        Qext_ += (n1 + n1 + 1.0) * (an_[n].real());
        Qsca_ += (n1 + n1 + 1.0) * (an_[n].real() * an_[n].real() +
                                    an_[n].imag() * an_[n].imag());
        Qpr_ += std::nan("");
        Qbktmp +=
            (FloatType)(n1 + n1 + 1.0) * (1.0 - 2.0 * (n1 % 2)) * (an_[n]);
        for (unsigned int t = 0; t < theta_.size(); t++) {
          calcPiTau(costheta[t], Pi, Tau);
          S1_[t] += calc_S1(n1, an_[n], static_cast<std::complex<FloatType>>(0),
                            Pi[n], Tau[n]);
          S2_[t] += calc_S2(n1, an_[n], static_cast<std::complex<FloatType>>(0),
                            Pi[n], Tau[n]);
        }
      }
      if (mode_type_ == Modes::kMagnetic || mode_type_ == Modes::kAll) {
        Qext_ += (n1 + n1 + 1.0) * (bn_[n].real());
        Qsca_ += (n1 + n1 + 1.0) * (bn_[n].real() * bn_[n].real() +
                                    bn_[n].imag() * bn_[n].imag());
        Qpr_ += std::nan("");
        Qbktmp +=
            (FloatType)(n1 + n1 + 1.0) * (1.0 - 2.0 * (n1 % 2)) * (bn_[n]);
        for (unsigned int t = 0; t < theta_.size(); t++) {
          calcPiTau(costheta[t], Pi, Tau);
          S1_[t] += calc_S1(n1, static_cast<std::complex<FloatType>>(0), bn_[n],
                            Pi[n], Tau[n]);
          S2_[t] += calc_S2(n1, static_cast<std::complex<FloatType>>(0), bn_[n],
                            Pi[n], Tau[n]);
        }
      }
    }
  }
  FloatType x2 = pow2(x.back());
  Qext_ = 2.0 * (Qext_) / x2;                  // Equation (27)
  Qsca_ = 2.0 * (Qsca_) / x2;                  // Equation (28)
  Qpr_ = Qext_ - 4.0 * (Qpr_) / x2;            // Equation (29)
  Qabs_ = Qext_ - Qsca_;                       // Equation (30)
  albedo_ = Qsca_ / Qext_;                     // Equation (31)
  asymmetry_factor_ = (Qext_ - Qpr_) / Qsca_;  // Equation (32)
  Qbk_ = (Qbktmp.real() * Qbktmp.real() + Qbktmp.imag() * Qbktmp.imag()) /
         x2;  // Equation (33)

  isMieCalculated_ = true;
}

}  // end of namespace nmie
#endif  // SRC_NMIE_BASIC_HPP_
