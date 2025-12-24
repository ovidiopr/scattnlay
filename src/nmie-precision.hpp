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
#endif

#ifdef MULTI_PRECISION
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/number.hpp>
#endif  // MULTI_PRECISION

namespace nmie {

#ifdef WITH_HWY
namespace hn = hwy::HWY_NAMESPACE;

template <typename T_base>
struct HighwayEngine {
  // Highway handles vectors of T_base (double or float)
  using T = T_base;
  using D = hn::ScalableTag<T>;
  using V = hn::Vec<D>;

  // Represents a batch of complex numbers (Real vectors, Imag vectors)
  struct ComplexV { V re; V im; };

  static inline V sin(V v) { 
    // Note: Highway provides math functions in hwy/contrib/math/math-inl.h
    // For smoke test, we use basic arithmetic
    return v; // Placeholder for actual vectorized sin
  }
  static inline V cos(V v) { (void)v; return hn::Undefined(D()); }
  static inline V exp(V v) { (void)v; return hn::Undefined(D()); }
  static inline V tan(V v) { (void)v; return hn::Undefined(D()); }
  
  static inline V add(V a, V b) {
    return hn::Add(a, b);
  }

  static inline V sub(V a, V b) {
    return hn::Sub(a, b);
  }

  static inline V mul(V a, V b) {
    return hn::Mul(a, b);
  }

  static inline V div(V a, V b) {
    return hn::Div(a, b);
  }

  static inline V set(T val) { return hn::Set(D(), val); }
  
  static inline V sign(V v) {
    auto zero = hn::Zero(D());
    auto one = hn::Set(D(), 1.0);
    auto minus_one = hn::Set(D(), -1.0);
    auto mask = hn::Gt(v, zero);
    return hn::IfThenElse(mask, one, minus_one);
  }

  // Complex arithmetic
  static inline ComplexV add(ComplexV a, ComplexV b) {
    return {hn::Add(a.re, b.re), hn::Add(a.im, b.im)};
  }

  static inline ComplexV sub(ComplexV a, ComplexV b) {
    return {hn::Sub(a.re, b.re), hn::Sub(a.im, b.im)};
  }

  static inline ComplexV mul(ComplexV a, ComplexV b) {
    // (a.re*b.re - a.im*b.im) + i(a.re*b.im + a.im*b.re)
    return {hn::Sub(hn::Mul(a.re, b.re), hn::Mul(a.im, b.im)),
            hn::Add(hn::Mul(a.re, b.im), hn::Mul(a.im, b.re))};
  }
  
  static inline ComplexV div(ComplexV a, ComplexV b) {
    V mag_sq = hn::Add(hn::Mul(b.re, b.re), hn::Mul(b.im, b.im));
    V re = hn::Div(hn::Add(hn::Mul(a.re, b.re), hn::Mul(a.im, b.im)), mag_sq);
    V im = hn::Div(hn::Sub(hn::Mul(a.im, b.re), hn::Mul(a.re, b.im)), mag_sq);
    return {re, im};
  }

  static inline V get_real(ComplexV z) { return z.re; }
  static inline V get_imag(ComplexV z) { return z.im; }
  static inline ComplexV make_complex(V re, V im) { return {re, im}; }
};
#endif

#ifdef MULTI_PRECISION
namespace nmm = boost::multiprecision;
typedef nmm::number<nmm::cpp_bin_float<MULTI_PRECISION> > FloatType;
#else
namespace nmm = std;
typedef double FloatType;
// typedef float FloatType;
#endif  // MULTI_PRECISION

struct ScalarEngine {
  using T = FloatType;
  static inline T sin(T v) { return nmm::sin(v); }
  static inline T cos(T v) { return nmm::cos(v); }
  static inline T exp(T v) { return nmm::exp(v); }
  static inline T sqrt(T v) { return nmm::sqrt(v); }
  static inline T abs(T v) { return nmm::abs(v); }
  static inline T log(T v) { return nmm::log(v); }
  static inline T tan(T v) { return nmm::tan(v); }
  static inline T ceil(T v) { return nmm::ceil(v); }
  static inline T max(T a, T b) { return std::max(a, b); }
  static inline T pow(T b, T e) { return nmm::pow(b, e); }

  static inline T add(T a, T b) { return a + b; }
  static inline T sub(T a, T b) { return a - b; }
  static inline T mul(T a, T b) { return a * b; }
  static inline T div(T a, T b) { return a / b; }

  static inline T set(T val) { return val; }
  static inline T sign(T v) { return (v > 0) ? 1 : -1; }

  static inline std::complex<T> add(std::complex<T> a, std::complex<T> b) { return a + b; }
  static inline std::complex<T> sub(std::complex<T> a, std::complex<T> b) { return a - b; }
  static inline std::complex<T> mul(std::complex<T> a, std::complex<T> b) { return a * b; }
  static inline std::complex<T> div(std::complex<T> a, std::complex<T> b) { return a / b; }

  static inline T get_real(std::complex<T> z) { return z.real(); }
  static inline T get_imag(std::complex<T> z) { return z.imag(); }
  static inline std::complex<T> make_complex(T re, T im) { return std::complex<T>(re, im); }
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
