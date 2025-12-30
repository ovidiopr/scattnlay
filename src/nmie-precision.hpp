#ifndef SRC_NMIE_PRECISION_H_
#define SRC_NMIE_PRECISION_H_
//**********************************************************************************//
//    Copyright (C) 2009-2018  Ovidio Pena <ovidio@bytesfall.com> // Copyright
//    (C) 2013-2018  Konstantin Ladutenko <kostyfisik@gmail.com>          //
//                                                                                  //
//    This file is part of scattnlay //
//                                                                                  //
//    This program is free software: you can redistribute it and/or modify // it
//    under the terms of the GNU General Public License as published by // the
//    Free Software Foundation, either version 3 of the License, or // (at your
//    option) any later version.                                           //
//                                                                                  //
//    This program is distributed in the hope that it will be useful, // but
//    WITHOUT ANY WARRANTY; without even the implied warranty of //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the // GNU
//    General Public License for more details. //
//                                                                                  //
//    The only additional remark is that we expect that all publications //
//    describing work using this software, or all commercial products // using
//    it, cite at least one of the following references:                      //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by //
//        a multilayered sphere," Computer Physics Communications, // vol. 180,
//        Nov. 2009, pp. 2348-2354.                                       //
//    [2] K. Ladutenko, U. Pal, A. Rivera, and O. Pena-Rodriguez, "Mie //
//        calculation of electromagnetic near-field for a multilayered //
//        sphere," Computer Physics Communications, vol. 214, May 2017, // pp.
//        225-230. //
//                                                                                  //
//    You should have received a copy of the GNU General Public License // along
//    with this program.  If not, see <http://www.gnu.org/licenses/>.         //
//**********************************************************************************//

//**********************************************************************************//
// This class implements the algorithm for a multilayered sphere described by:
// //
//    [1] W. Yang, "Improved recursive algorithm for light scattering by a //
//        multilayered sphere,” Applied Optics, vol. 42, Mar. 2003, pp.
//        1710-1720.  //
//                                                                                  //
// You can find the description of all the used equations in: //
//    [2] O. Pena and U. Pal, "Scattering of electromagnetic radiation by //
//        a multilayered sphere," Computer Physics Communications, // vol. 180,
//        Nov. 2009, pp. 2348-2354.                                       //
//                                                                                  //
// Hereinafter all equations numbers refer to [2] //
//**********************************************************************************//
#include <complex>
#include <vector>

#ifdef WITH_HWY
#include "hwy/highway.h"
#include "hwy/contrib/math/math-inl.h"
#endif

#ifdef MULTI_PRECISION
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/number.hpp>
#endif  // MULTI_PRECISION

namespace nmie {



#ifdef WITH_HWY
namespace hn = hwy::HWY_NAMESPACE;
#include "highway-engine.hpp"
#endif

#ifdef MULTI_PRECISION
namespace nmm = boost::multiprecision;
typedef nmm::number<nmm::cpp_bin_float<MULTI_PRECISION> > FloatType;
#else
namespace nmm = std;
typedef double FloatType;
// typedef float FloatType;
#endif  // MULTI_PRECISION

template <typename T>
struct ScalarEngine {
  using RealV = T;
  using ComplexV = std::complex<T>;
  using MaskV = bool; // For scalar, masks are just booleans

  // Traits
  static constexpr size_t Lanes() { return 1; }

  // Memory Operations
  static inline RealV set(T val) { return val; }
  static inline void store_interleaved(ComplexV z, std::complex<T>* ptr) { *ptr = z; }
  static inline ComplexV load_interleaved(const std::complex<T>* ptr) { return *ptr; }

  static inline void store(ComplexV z, std::complex<T>* ptr) { *ptr = z; }
  static inline ComplexV load(const std::complex<T>* ptr) { return *ptr; }

  static inline void store(RealV v, T* ptr) { *ptr = v; }
  static inline RealV load(const T* ptr) { return *ptr; }

  // Math Primitives (Using ADL to handle std:: vs nmm::)
  static inline RealV abs(RealV v) { using std::abs; using nmm::abs; return abs(v); }
  static inline RealV sin(RealV v) { using std::sin; using nmm::sin; return sin(v); }
  static inline RealV cos(RealV v) { using std::cos; using nmm::cos; return cos(v); }
  static inline RealV exp(RealV v) { using std::exp; using nmm::exp; return exp(v); }
  static inline RealV sqrt(RealV v) { using std::sqrt; using nmm::sqrt; return sqrt(v); }
  static inline RealV log(RealV v) { using std::log; using nmm::log; return log(v); }
  static inline RealV tan(RealV v) { using std::tan; using nmm::tan; return tan(v); }
  static inline RealV ceil(RealV v) { using std::ceil; using nmm::ceil; return ceil(v); }
  static inline RealV floor(RealV v) { using std::floor; using nmm::floor; return floor(v); }
  static inline RealV max(RealV a, RealV b) { using std::max; return max(a, b); }
  static inline RealV pow(RealV b, RealV e) { using std::pow; using nmm::pow; return pow(b, e); }

  static inline RealV add(RealV a, RealV b) { return a + b; }
  static inline RealV sub(RealV a, RealV b) { return a - b; }
  static inline RealV mul(RealV a, RealV b) { return a * b; }
  static inline RealV div(RealV a, RealV b) { return a / b; }

  static inline RealV sign(RealV v) { return (v > 0) ? 1 : -1; }

  static inline ComplexV add(ComplexV a, ComplexV b) { return a + b; }
  static inline ComplexV sub(ComplexV a, ComplexV b) { return a - b; }
  static inline ComplexV mul(ComplexV a, ComplexV b) { return a * b; }
  static inline ComplexV div(ComplexV a, ComplexV b) { return a / b; }

  static inline RealV get_real(ComplexV z) { return z.real(); }
  static inline RealV get_imag(ComplexV z) { return z.imag(); }
  static inline ComplexV make_complex(RealV re, RealV im) { return ComplexV(re, im); }

  static inline ComplexV log(ComplexV z) {
    RealV re = z.real();
    RealV im = z.imag();
    using std::sqrt; using nmm::sqrt;
    using std::atan2; using nmm::atan2;
    using std::log; using nmm::log;
    RealV r = sqrt(re * re + im * im);
    RealV phi = atan2(im, re);
    return ComplexV(log(r), phi);
  }

  static inline ComplexV sqrt(ComplexV z) {
    RealV re = z.real();
    RealV im = z.imag();
    using std::sqrt; using nmm::sqrt;
    RealV r = sqrt(re * re + im * im);
    RealV u = sqrt((r + re) * T(0.5));
    RealV v = sqrt((r - re) * T(0.5));
    if (im < 0) v = -v;
    return ComplexV(u, v);
  }

  static inline RealV reduce_max(RealV v) { return v; }
};

template <class T>
T sin_t(T v) {
  if (std::is_same<T, double>::value)
    return static_cast<T>(std::sin(static_cast<double>(v)));
  return static_cast<T>(nmm::sin(static_cast<FloatType>(v)));
}
template <class T>
T cos_t(T v) {
  if (std::is_same<T, double>::value)
    return static_cast<T>(std::cos(static_cast<double>(v)));
  return static_cast<T>(nmm::cos(static_cast<FloatType>(v)));
}
template <class T>
T sqrt_t(T v) {
  if (std::is_same<T, double>::value)
    return static_cast<T>(std::sqrt(static_cast<double>(v)));
  return static_cast<T>(nmm::sqrt(static_cast<FloatType>(v)));
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
std::vector<std::complex<ToFloatType> > ConvertComplexVector(
    std::vector<std::complex<FromFloatType> > x) {
  std::vector<std::complex<ToFloatType> > new_x;
  for (auto element : x) {
    new_x.push_back(
        std::complex<ToFloatType>(static_cast<ToFloatType>(element.real()),
                                  static_cast<ToFloatType>(element.imag())));
  }
  return new_x;
}

template <typename ToFloatType, typename FromFloatType>
std::vector<std::vector<std::complex<ToFloatType> > >
ConvertComplexVectorVector(
    std::vector<std::vector<std::complex<FromFloatType> > > x) {
  std::vector<std::vector<std::complex<ToFloatType> > > new_x;
  std::vector<std::complex<ToFloatType> > new_y;
  for (auto y : x) {
    new_y.clear();
    for (auto element : y) {
      new_y.push_back(
          std::complex<ToFloatType>(static_cast<ToFloatType>(element.real()),
                                    static_cast<ToFloatType>(element.imag())));
    }
    new_x.push_back(new_y);
  }
  return new_x;
}

}  // end of namespace nmie
#endif  // SRC_NMIE_PRECISION_H_
