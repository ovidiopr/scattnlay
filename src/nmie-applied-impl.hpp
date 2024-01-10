#ifndef SRC_NMIE_APPLIED_IMPL_HPP_
#define SRC_NMIE_APPLIED_IMPL_HPP_
//**********************************************************************************
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
//
//    @brief  Wrapper class around nMie function for ease of use
//
//**********************************************************************************
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>

namespace nmie {
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::GetFailed() {
  FloatType faild_x = 9.42477796076938;
  // FloatType faild_x = 9.42477796076937;
  std::complex<FloatType> z(faild_x, 0.0);
  std::vector<int> nmax_local_array = {20, 100, 500, 2500};
  for (auto nmax_local : nmax_local_array) {
    std::vector<std::complex<FloatType> > D1_failed(nmax_local + 1);
    // Downward recurrence for D1 - equations (16a) and (16b)
    D1_failed[nmax_local] = std::complex<FloatType>(0.0, 0.0);
    const std::complex<FloatType> zinv = std::complex<FloatType>(1.0, 0.0) / z;
    for (int n = nmax_local; n > 0; n--) {
      D1_failed[n - 1] =
          FloatType(n) * zinv - 1.0 / (D1_failed[n] + FloatType(n) * zinv);
    }
    printf("Faild D1[0] from reccurence (z = %16.14f, nmax = %d): %g\n",
           faild_x, nmax_local, D1_failed[0].real());
  }
  printf("Faild D1[0] from continued fraction (z = %16.14f): %g\n", faild_x,
         calcD1confra(0, z).real());
  // D1[nmax_] = calcD1confra(nmax_, z);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::AddTargetLayer(
    FloatType width,
    std::complex<FloatType> layer_index) {
  this->MarkUncalculated();
  if (width <= 0)
    throw std::invalid_argument("Layer width should be positive!");
  target_width_.push_back(width);
  target_index_.push_back(layer_index);
}  // end of void  MultiLayerMieApplied<FloatType>::AddTargetLayer(...)
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::AddTargetLayerReIm(
    FloatType width,
    FloatType re_layer_index,
    FloatType im_layer_index) {
  this->MarkUncalculated();
  if (width <= 0)
    throw std::invalid_argument("Layer width should be positive!");
  target_width_.push_back(width);
  target_index_.push_back(
      std::complex<FloatType>(re_layer_index, im_layer_index));
}  // end of void  MultiLayerMieApplied<FloatType>::AddTargetLayer(...)
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::SetTargetPEC(FloatType radius) {
  this->MarkUncalculated();
  if (target_width_.size() != 0 || target_index_.size() != 0)
    throw std::invalid_argument(
        "Error! Define PEC target radius before any other layers!");
  // Add layer of any index...
  AddTargetLayer(radius, std::complex<FloatType>(0.0, 0.0));
  // ... and mark it as PEC
  this->SetPECLayer(0);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::SetCoatingIndex(
    std::vector<std::complex<FloatType> > index) {
  this->MarkUncalculated();
  coating_index_.clear();
  for (auto value : index)
    coating_index_.push_back(value);
}  // end of void
   // MultiLayerMieApplied<FloatType>::SetCoatingIndex(std::vector<complex>
   // index);
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::SetCoatingWidth(
    std::vector<FloatType> width) {
  this->MarkUncalculated();
  coating_width_.clear();
  for (auto w : width)
    if (w <= 0)
      throw std::invalid_argument("Coating width should be positive!");
    else
      coating_width_.push_back(w);
}
// end of void MultiLayerMieApplied<FloatType>::SetCoatingWidth(...);
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::SetWidthSP(
    const std::vector<FloatType>& size_parameter) {
  this->MarkUncalculated();
  this->size_param_.clear();
  FloatType prev_size_parameter = 0.0;
  for (auto layer_size_parameter : size_parameter) {
    if (layer_size_parameter <= 0.0)
      throw std::invalid_argument("Size parameter should be positive!");
    if (prev_size_parameter > layer_size_parameter)
      throw std::invalid_argument(
          "Size parameter for next layer should be larger than the previous "
          "one!");
    prev_size_parameter = layer_size_parameter;
    this->size_param_.push_back(layer_size_parameter);
  }
}
// end of void MultiLayerMieApplied<FloatType>::SetWidthSP(...);
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::SetIndexSP(
    const std::vector<std::complex<FloatType> >& index) {
  this->MarkUncalculated();
  // refractive_index_.clear();
  this->refractive_index_ = index;
  // for (auto value : index) refractive_index_.push_back(value);
}  // end of void MultiLayerMieApplied<FloatType>::SetIndexSP(...);
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::SetFieldPointsSP(
    const std::vector<std::vector<FloatType> >& coords_sp) {
  if (coords_sp.size() != 3)
    throw std::invalid_argument(
        "Error! Wrong dimension of field monitor points!");
  if (coords_sp[0].size() != coords_sp[1].size() ||
      coords_sp[0].size() != coords_sp[2].size())
    throw std::invalid_argument(
        "Error! Missing coordinates for field monitor points!");
  this->coords_ = coords_sp;
  // for (int i = 0; i < coords_sp_[0].size(); ++i) {
  //   printf("%g, %g, %g\n", coords_sp_[0][i], coords_sp_[1][i],
  //   coords_sp_[2][i]);
  // }
}  // end of void MultiLayerMieApplied<FloatType>::SetFieldPointsSP(...)
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::GenerateSizeParameter() {
  this->MarkUncalculated();
  this->size_param_.clear();
  FloatType radius = 0.0;
  for (auto width : target_width_) {
    radius += width;
    this->size_param_.push_back(2 * nmie::PI_ * radius / wavelength_);
  }
  for (auto width : coating_width_) {
    radius += width;
    this->size_param_.push_back(2 * nmie::PI_ * radius / wavelength_);
  }
  this->total_radius_ = radius;
}  // end of void MultiLayerMieApplied<FloatType>::GenerateSizeParameter();
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::GenerateIndex() {
  this->MarkUncalculated();
  this->refractive_index_.clear();
  for (auto index : this->target_index_)
    this->refractive_index_.push_back(index);
  for (auto index : this->coating_index_)
    this->refractive_index_.push_back(index);
}  // end of void MultiLayerMieApplied<FloatType>::GenerateIndex();
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
FloatType MultiLayerMieApplied<FloatType>::GetTotalRadius() {
  if (!this->isMieCalculated())
    GenerateSizeParameter();
  return this->total_radius_;
}  // end of FloatType MultiLayerMieApplied<FloatType>::GetTotalRadius();
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
std::vector<std::vector<FloatType> >
MultiLayerMieApplied<FloatType>::GetSpectra(FloatType from_WL,
                                            FloatType to_WL,
                                            int samples) {
  if (!this->isMieCalculated())
    throw std::invalid_argument(
        "You should run calculations before result request!");
  std::vector<std::vector<FloatType> > spectra;
  FloatType step_WL = (to_WL - from_WL) / static_cast<FloatType>(samples);
  FloatType wavelength_backup = wavelength_;
  long fails = 0;
  for (FloatType WL = from_WL; WL < to_WL; WL += step_WL) {
    wavelength_ = WL;
    try {
      RunMieCalculation();
    } catch (const std::invalid_argument& ia) {
      fails++;
      continue;
    }
    // printf("%3.1f ",WL);
    spectra.push_back(
        std::vector<FloatType>({wavelength_, this->GetQext(), this->GetQsca(),
                                this->GetQabs(), this->GetQbk()}));
  }  // end of for each WL in spectra
  printf("Spectrum has %li fails\n", fails);
  wavelength_ = wavelength_backup;
  return spectra;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::ClearTarget() {
  this->MarkUncalculated();
  this->target_width_.clear();
  this->target_index_.clear();
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::ClearCoating() {
  this->MarkUncalculated();
  this->coating_width_.clear();
  this->coating_index_.clear();
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::ClearLayers() {
  this->MarkUncalculated();
  this->ClearTarget();
  this->ClearCoating();
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::ClearAllDesign() {
  this->MarkUncalculated();
  this->ClearLayers();
  this->size_param_.clear();
  this->refractive_index_.clear();
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
//                         Computational core
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::ConvertToSP() {
  this->MarkUncalculated();
  if (target_width_.size() + coating_width_.size() == 0)
    return;  // Nothing to convert, we suppose that SP was set directly
  GenerateSizeParameter();
  GenerateIndex();
  if (this->size_param_.size() != this->refractive_index_.size())
    throw std::invalid_argument(
        "Ivalid conversion of width to size parameter units!/n");
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::RunMieCalculation() {
  ConvertToSP();
  this->MultiLayerMie<FloatType>::RunMieCalculation();
}

template <typename FloatType>
void MultiLayerMieApplied<FloatType>::RunFieldCalculationPolar(
    const int outer_arc_points,
    const int radius_points,
    const double from_Rho,
    const double to_Rho,
    const double from_Theta,
    const double to_Theta,
    const double from_Phi,
    const double to_Phi,
    const int isIgnoreAvailableNmax) {
  ConvertToSP();  // Converts to size parameter units only the particle design,
  // so we need to convert input parameters too...
  const FloatType a = 2 * nmie::PI_ / wavelength_;
  this->MultiLayerMie<FloatType>::RunFieldCalculationPolar(
      outer_arc_points, radius_points, a * from_Rho, a * to_Rho, from_Theta,
      to_Theta, from_Phi, to_Phi, isIgnoreAvailableNmax == 0 ? false : true);
}

template <typename FloatType>
void MultiLayerMieApplied<FloatType>::RunFieldCalculationCartesian(
    const int first_side_points,
    const int second_side_points,
    const double relative_side_length,
    const int plane_selected,
    const double at_x,
    const double at_y,
    const double at_z,
    const int isIgnoreAvailableNmax) {
  //  std::cout<<'test'<<std::endl;
  ConvertToSP();
  this->MultiLayerMie<FloatType>::RunFieldCalculationCartesian(
      first_side_points, second_side_points, relative_side_length,
      plane_selected, at_x, at_y, at_z,
      isIgnoreAvailableNmax == 0 ? false : true);
}

// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
void MultiLayerMieApplied<FloatType>::GetExpanCoeffs(
    std::vector<std::vector<std::complex<FloatType> > >& aln,
    std::vector<std::vector<std::complex<FloatType> > >& bln,
    std::vector<std::vector<std::complex<FloatType> > >& cln,
    std::vector<std::vector<std::complex<FloatType> > >& dln) {
  ConvertToSP();  // Case of call before running full Mie calculation.
  // Calculate scattering coefficients an_ and bn_
  this->calcScattCoeffs();
  // Calculate expansion coefficients aln_,  bln_, cln_, and dln_
  this->calcExpanCoeffs();
  aln = this->aln_;
  bln = this->bln_;
  cln = this->cln_;
  dln = this->dln_;

}  // end of void MultiLayerMieApplied<FloatType>::GetExpanCoeffs( ...)
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //

}  // end of namespace nmie
#endif  // SRC_NMIE_APPLIED_IMPL_HPP_
