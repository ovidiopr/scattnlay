#ifndef SRC_NMIE_WEB_HPP_
#define SRC_NMIE_WEB_HPP_
//**********************************************************************************
//    Copyright (C) 2009-2024  Ovidio Pena <ovidio@bytesfall.com>
//    Copyright (C) 2013-2024  Konstantin Ladutenko <kostyfisik@gmail.com>
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
//**********************************************************************************

#include <array>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "nmie-applied-impl.hpp"

#include <emscripten/bind.h>
#include <emscripten/val.h>

namespace nmie {

template <typename FloatType = double>
class MultiLayerMieWeb : public MultiLayerMieApplied<FloatType> {
  // Will throw for any error!
 public:
  emscripten::val GetFieldEabs();
  void ClearTarget() { this->MultiLayerMieApplied<FloatType>::ClearTarget(); };
  void SetWavelength(FloatType wavelength) {
    this->MultiLayerMieApplied<FloatType>::SetWavelength(wavelength);
  };
  void SetModeNmaxAndType(int mode_n, int mode_type) {
    this->MultiLayerMieApplied<FloatType>::SetModeNmaxAndType(mode_n,
                                                              mode_type);
  };
  void AddTargetLayerReIm(FloatType width,
                          FloatType re_layer_index,
                          FloatType im_layer_index) {
    this->MultiLayerMieApplied<FloatType>::AddTargetLayerReIm(
        width, re_layer_index, im_layer_index);
  };
  void RunMieCalculation() {
    this->MultiLayerMieApplied<FloatType>::RunMieCalculation();
  };
  void RunFieldCalculationPolar(const int outer_arc_points,
                                const int radius_points,
                                const double from_Rho,
                                const double to_Rho,
                                const double from_Theta,
                                const double to_Theta,
                                const double from_Phi,
                                const double to_Phi,
                                const int isIgnoreAvailableNmax) {
    this->MultiLayerMieApplied<FloatType>::RunFieldCalculationPolar(
        outer_arc_points, radius_points, from_Rho, to_Rho, from_Theta, to_Theta,
        from_Phi, to_Phi, isIgnoreAvailableNmax);
  };
  void RunFieldCalculationCartesian(const int first_side_points,
                                    const int second_side_points,
                                    const double relative_side_length,
                                    const int plane_selected,
                                    const double at_x,
                                    const double at_y,
                                    const double at_z,
                                    const int isIgnoreAvailableNmax) {
    this->MultiLayerMieApplied<FloatType>::RunFieldCalculationCartesian(
        first_side_points, second_side_points, relative_side_length,
        plane_selected, at_x, at_y, at_z, isIgnoreAvailableNmax);
  };
  FloatType GetQsca() {
    return this->MultiLayerMieApplied<FloatType>::GetQsca();
  };
  FloatType GetQabs() {
    return this->MultiLayerMieApplied<FloatType>::GetQabs();
  };
  FloatType GetQext() {
    return this->MultiLayerMieApplied<FloatType>::GetQext();
  };

};  // end of class MultiLayerMieWeb

}  // end of namespace nmie
#endif  // SRC_NMIE_WEB_HPP
