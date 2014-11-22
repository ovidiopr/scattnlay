#ifndef SRC_NMIE_NMIE_WRAPPER_H_
#define SRC_NMIE_NMIE_WRAPPER_H_
///
/// @file   nmie-wrapper.h
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Tue Sep  3 00:40:47 2013
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
///
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::exit(EXIT_FAILURE); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

namespace nmie {
  class MultiLayerMie {
   public:
    void SetWavelength(double wavelength) {wavelength_ = wavelength;};
    void SetSizeParameter(std::vector<double> size_parameter);
    void SetIndex(std::vector< std::complex<double> > index);
    void RunMie(double *Qext, double *Qsca, double *Qabs, double *Qbk);
    std::vector< std::vector<double> >  GetSpectra(double from_WL, double to_WL,
                                                   int samples);
    /// Disabled functions or unused untill the end of C++ conversion
    double GetTotalRadius();
  private:
    /// Size parameter for each layer.
    std::vector<double> size_parameter_;
    /// Complex index value for each layer.
    std::vector< std::complex<double> > index_;
    const double PI=3.14159265358979323846;
    /// Disabled functions or unused untill the end of C++ conversion
    void GenerateSizeParameter();
    double wavelength_ = 1.0;
    double total_radius_ = 0.0;

  };  // end of class MultiLayerMie
}  // end of namespace nmie
#endif  // SRC_NMIE_NMIE_WRAPPER_H_
