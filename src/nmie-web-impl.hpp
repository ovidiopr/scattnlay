#ifndef SRC_NMIE_WEB_IMPL_HPP_
#define SRC_NMIE_WEB_IMPL_HPP_
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

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include "nmie-web.hpp"

namespace nmie {
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //

// from https://toughengineer.github.io/demo/dsp/fft-perf/
template <typename FloatType = double>
emscripten::val toJSFloat64Array(const std::vector<double>& v) {
  emscripten::val view{emscripten::typed_memory_view(
      v.size(), v.data())};  // make a view of local object

  auto result = emscripten::val::global("Float64Array")
                    .new_(v.size());  // make a JS typed array
  result.call<void>(
      "set",
      view);  // set typed array values "on the JS side" using the memory view

  return result;
}

// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
template <typename FloatType>
emscripten::val MultiLayerMieWeb<FloatType>::GetFieldEabs() {
  auto Eabs = this->MultiLayerMie<FloatType>::GetFieldEabs();
  return toJSFloat64Array(Eabs);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
}  // end of namespace nmie
#endif  // SRC_NMIE_WEB_IMPL_HPP_
