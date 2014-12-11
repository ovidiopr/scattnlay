///
/// @file   nmie-wrapper.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Tue Sep  3 00:38:27 2013
/// @copyright 2013 Ladutenko Konstantin
///
/// nmie-wrapper is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// nmie-wrapper is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with nmie-wrapper.  If not, see <http://www.gnu.org/licenses/>.
///
/// nmie-wrapper uses nmie.c from scattnlay by Ovidio Pena
/// <ovidio@bytesfall.com> as a linked library. He has an additional condition to 
/// his library:
//    The only additional condition is that we expect that all publications         //
//    describing  work using this software , or all commercial products             //
//    using it, cite the following reference:                                       //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
///
/// @brief  Wrapper class around nMie function for ease of use
///
#include "ucomplex.h"
#include "nmie-wrapper.h"
#include "nmie.h"
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>

namespace nmie {  
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
    void MultiLayerMie::SetSizeParameter(std::vector<double> size_parameter) {
    size_parameter_.clear();
    for (auto layer_size_parameter : size_parameter)
      if (layer_size_parameter <= 0.0)
        throw std::invalid_argument("Size parameter should be positive!");
      else size_parameter_.push_back(layer_size_parameter);
  }
  // end of void MultiLayerMie::SetSizeParameter(...);
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::SetIndex(std::vector< std::complex<double> > index) {
    index_.clear();
    for (auto value : index) index_.push_back(value);
  }  // end of void MultiLayerMie::SetIndex(...);
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// nMie starts layer numeration from 1 (no separation of target
  /// and coating). So the first elment (zero-indexed) is not used
  /// and has some unused value.
  /// Kostya, that's no longer the case. Now the layer numbering starts at zero.
  void MultiLayerMie::GenerateSizeParameter() {
    // size_parameter_.clear();
    // size_parameter_.push_back(0.0);
    // double radius = 0.0;
    // for (auto width : target_thickness_) {
    //   radius += width;
    //   size_parameter_.push_back(2*PI*radius / wavelength_);
    // }
    // for (auto width : coating_thickness_) {
    //   radius += width;
    //   size_parameter_.push_back(2*PI*radius / wavelength_);
    // }
    // total_radius_ = radius;
  }  // end of void MultiLayerMie::GenerateSizeParameter();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::RunMie(double *Qext_out, double *Qsca_out,
                             double *Qabs_out, double *Qbk_out) {
    if (size_parameter_.size() != index_.size())
      throw std::invalid_argument("Each size parameter should have only one index!");
    int L = static_cast<int>(size_parameter_.size()) - 1;
    double Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo;
    int nt = 0;
    double *Theta = (double *)malloc(nt * sizeof(double));
    complex *S1 = (complex *)malloc(nt * sizeof(complex *));
    complex *S2 = (complex *)malloc(nt * sizeof(complex *));
    double *x = &(size_parameter_.front());
    complex *m = &(index_.front());
    int terms = 0;
    terms = nMie(L, x, m, nt, Theta,
                 &Qext, &Qsca, &Qabs, &Qbk, &Qpr, &g, &Albedo,
                 S1,S2);
    free(Theta);
    free(S1);
    free(S2);
    if (terms == 0) {
      *Qext_out = Qfaild_;
      *Qsca_out = Qfaild_;
      *Qabs_out = Qfaild_;
      *Qbk_out = Qfaild_;
      throw std::invalid_argument("Failed to evaluate Q!");
    }
    *Qext_out = Qext;
    *Qsca_out = Qsca;
    *Qabs_out = Qabs;
    *Qbk_out = Qbk;
  }  // end of void MultiLayerMie::RunMie();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double MultiLayerMie::GetTotalRadius() {
    if (total_radius_ == 0) GenerateSizeParameter();
    return total_radius_;      
  }  // end of double MultiLayerMie::GetTotalRadius();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector< std::vector<double> >
  MultiLayerMie::GetSpectra(double from_WL, double to_WL, int samples) {
    std::vector< std::vector<double> > spectra;
    double step_WL = (to_WL - from_WL)/ static_cast<double>(samples);
    double wavelength_backup = wavelength_;
    long fails = 0;
    for (double WL = from_WL; WL < to_WL; WL += step_WL) {
      double Qext, Qsca, Qabs, Qbk;
      wavelength_ = WL;
      try {
        RunMie(&Qext, &Qsca, &Qabs, &Qbk);
      } catch( const std::invalid_argument& ia ) {
        fails++;
        continue;
      }
      //printf("%3.1f ",WL);
      spectra.push_back({wavelength_, Qext, Qsca, Qabs, Qbk});
    }  // end of for each WL in spectra
    printf("fails %li\n",fails);
    wavelength_ = wavelength_backup;
    return spectra;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
///MultiLayerMie::
}  // end of namespace nmie
