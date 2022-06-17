#ifndef SRC_PB11_HELPERS_HPP_
#define SRC_PB11_HELPERS_HPP_
//******************************************************************************
//    Copyright (C) 2009-2022  Ovidio Pena <ovidio@bytesfall.com>
//    Copyright (C) 2013-2022  Konstantin Ladutenko <kostyfisik@gmail.com>
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
//******************************************************************************

#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
//******************************************************************************
template <typename T>
std::vector<T> Py2Vector(const py::array_t<T>& py_x) {
  std::vector<T> c_x(py_x.size());
  std::memcpy(c_x.data(), py_x.data(), py_x.size() * sizeof(T));
  return c_x;
}

//******************************************************************************
// https://github.com/pybind/pybind11/issues/1042#issuecomment-508582847
// template <typename Sequence>
// inline py::array_t<typename Sequence::value_type> Vector2Py(Sequence&& seq) {
//  // Move entire object to heap (Ensure is moveable!). Memory handled via
//  Python capsule Sequence* seq_ptr = new Sequence(std::move(seq)); auto
//  capsule = py::capsule(seq_ptr, [](void* p) { delete
//  reinterpret_cast<Sequence*>(p); }); return py::array(seq_ptr->size(),  //
//  shape of array
//                   seq_ptr->data(),  // c-style contiguous strides for
//                   Sequence capsule           // numpy array references this
//                   parent
//  );
//}

//******************************************************************************
template <typename outputType>
inline py::array_t<outputType> Vector2Py(const std::vector<outputType>& seq) {
  return py::array(seq.size(), seq.data());
}

//******************************************************************************
template <typename inputType = double, typename outputType = double>
py::array_t<std::complex<outputType>> VectorComplex2Py(
    const std::vector<std::complex<inputType>>& cf_x) {
  auto c_x = nmie::ConvertComplexVector<outputType, inputType>(cf_x);
  auto py_x = py::array_t<std::complex<outputType>>(c_x.size());
  auto py_x_buffer = py_x.request();
  auto* py_x_ptr = (std::complex<outputType>*)py_x_buffer.ptr;
  std::memcpy(py_x_ptr, c_x.data(),
              c_x.size() * sizeof(std::complex<outputType>));
  return py_x;
}

//******************************************************************************
// https://stackoverflow.com/questions/17294629/merging-flattening-sub-vectors-into-a-single-vector-c-converting-2d-to-1d
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v) {
  std::size_t total_size = 0;
  for (const auto& sub : v)
    total_size += sub.size();  // I wish there was a transform_accumulate
  std::vector<T> result;
  result.reserve(total_size);
  for (const auto& sub : v)
    result.insert(result.end(), sub.begin(), sub.end());
  return result;
}

//******************************************************************************
template <typename T>
py::array Vector2DComplex2Py(const std::vector<std::vector<T>>& x) {
  size_t dim1 = x.size();
  size_t dim2 = x[0].size();
  auto result = flatten(x);
  // https://github.com/tdegeus/pybind11_examples/blob/master/04_numpy-2D_cpp-vector/example.cpp
  size_t ndim = 2;
  std::vector<size_t> shape = {dim1, dim2};
  std::vector<size_t> strides = {sizeof(T) * dim2, sizeof(T)};

  // return 2-D NumPy array
  return py::array(py::buffer_info(
      result.data(),                      /* data as contiguous array  */
      sizeof(T),                          /* size of one scalar        */
      py::format_descriptor<T>::format(), /* data type                 */
      ndim,                               /* number of dimensions      */
      shape,                              /* shape of the matrix       */
      strides                             /* strides for each axis     */
      ));
}

#endif