#ifndef SRC_NMIE_PRECISION_H_
#define SRC_NMIE_PRECISION_H_
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
//                                                                                  //
// Hereinafter all equations numbers refer to [2]                                   //
//**********************************************************************************//
#ifdef MULTI_PRECISION
#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#endif  // MULTI_PRECISION

namespace nmie {
  #ifdef MULTI_PRECISION
  namespace nmm = boost::multiprecision;
  typedef nmm::number<nmm::cpp_bin_float<MULTI_PRECISION> > FloatType;
  #else
  namespace nmm = std;
  typedef double FloatType;
  //typedef float FloatType;
  #endif  // MULTI_PRECISION

  template<class T> T sin_t(T v) {
    if (std::is_same<T, double>::value) return static_cast<T>(std::sin(static_cast<double>(v)));
    return static_cast<T>(nmm::sin(static_cast<FloatType >(v)));
  }
  template<class T> T cos_t(T v) {
    if (std::is_same<T, double>::value) return static_cast<T>(std::cos(static_cast<double>(v)));
    return static_cast<T>(nmm::cos(static_cast<FloatType >(v)));
  }


template <typename ToFloatType, typename FromFloatType>
  std::vector<ToFloatType> ConvertVector(const std::vector<FromFloatType> x) {
    std::vector<ToFloatType> new_x;
    for (auto element : x) {
      new_x.push_back(static_cast<ToFloatType>(element));
    }
    return new_x;
  }

  template <typename ToFloatType, typename FromFloatType>
  std::complex<ToFloatType> ConvertComplex(std::complex<FromFloatType> z) {
    return std::complex<ToFloatType>(static_cast<ToFloatType>(z.real()),
                                     static_cast<ToFloatType>(z.imag()));
  }

 template <typename ToFloatType, typename FromFloatType>
  std::vector<std::complex<ToFloatType> > ConvertComplexVector(std::vector<std::complex<FromFloatType> > x) {
    std::vector<std::complex<ToFloatType> > new_x;
    for (auto element : x) {
      new_x.push_back(std::complex<ToFloatType>(static_cast<ToFloatType>(element.real()),
                                                static_cast<ToFloatType>(element.imag()) ) );
    }
    return new_x;
  }

  template <typename ToFloatType, typename FromFloatType>
  std::vector<std::vector<std::complex<ToFloatType> > > ConvertComplexVectorVector(std::vector<std::vector<std::complex<FromFloatType> > > x) {
    std::vector<std::vector<std::complex<ToFloatType> > > new_x;
    std::vector<std::complex<ToFloatType> >  new_y;
    for (auto y : x) {
      new_y.clear();
      for (auto element : y) {
        new_y.push_back(std::complex<ToFloatType>(static_cast<ToFloatType>(element.real()),
                                                  static_cast<ToFloatType>(element.imag()) ) );
      }
      new_x.push_back(new_y);
    }
    return new_x;
  }

}  // end of namespace nmie
#endif  // SRC_NMIE_PRECISION_H_
