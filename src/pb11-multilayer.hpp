#ifndef SRC_PB11_MULTILAYER_HPP_
#define SRC_PB11_MULTILAYER_HPP_
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

#include "pb11-helpers.hpp"

#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace nmie {
//******************************************************************************
// Python class interface declaration
//******************************************************************************
template <typename FloatType>
class PyMultiLayerMie : public MultiLayerMie<FloatType> {
 public:
  template <typename outputType>
  py::array_t<std::complex<outputType>> GetS1();
  template <typename outputType>
  py::array_t<std::complex<outputType>> GetS2();
  template <typename outputType>
  py::array_t<std::complex<outputType>> GetAn();
  template <typename outputType>
  py::array_t<std::complex<outputType>> GetBn();
  template <typename outputType>
  py::array GetLayerAn();
  template <typename outputType>
  py::array GetLayerBn();
  template <typename outputType>
  py::array GetLayerCn();
  template <typename outputType>
  py::array GetLayerDn();
  template <typename outputType>
  py::array GetFieldE();
  template <typename outputType>
  py::array GetFieldH();
  template <typename outputType>
  py::array_t<outputType> GetFieldEabs();
  template <typename outputType>
  py::array_t<outputType> GetFieldHabs();

  template <typename inputType>
  void SetLayersSize(
      const py::array_t<inputType, py::array::c_style | py::array::forcecast>&
          py_layer_size);
  template <typename inputType>
  void SetLayersIndex(
      const py::array_t<std::complex<inputType>,
                        py::array::c_style | py::array::forcecast>& py_index);
  template <typename inputType>
  void SetAngles(
      const py::array_t<inputType, py::array::c_style | py::array::forcecast>&
          py_angles);
  void SetFieldCoords(
      const py::array_t<double, py::array::c_style | py::array::forcecast>&
          py_Xp,
      const py::array_t<double, py::array::c_style | py::array::forcecast>&
          py_Yp,
      const py::array_t<double, py::array::c_style | py::array::forcecast>&
          py_Zp);
};

//******************************************************************************
// Python class interface implementation
//******************************************************************************
template <typename FloatType>
template <typename inputType>
void PyMultiLayerMie<FloatType>::SetLayersSize(
    const py::array_t<inputType, py::array::c_style | py::array::forcecast>&
        py_layer_size) {
  auto layer_size_dp = Py2Vector<inputType>(py_layer_size);
  this->MultiLayerMie<FloatType>::SetLayersSize(
      ConvertVector<FloatType>(layer_size_dp));
}

//******************************************************************************
template <typename FloatType>
template <typename inputType>
void PyMultiLayerMie<FloatType>::SetLayersIndex(
    const py::array_t<std::complex<inputType>,
                      py::array::c_style | py::array::forcecast>& py_index) {
  auto index_dp = Py2Vector<std::complex<inputType>>(py_index);
  this->MultiLayerMie<FloatType>::SetLayersIndex(
      ConvertComplexVector<FloatType>(index_dp));
}

//******************************************************************************
template <typename FloatType>
template <typename inputType>
void PyMultiLayerMie<FloatType>::SetAngles(
    const py::array_t<inputType, py::array::c_style | py::array::forcecast>&
        py_angles) {
  auto angles_dp = Py2Vector<inputType>(py_angles);
  this->MultiLayerMie<FloatType>::SetAngles(
      ConvertVector<FloatType>(angles_dp));
}

//******************************************************************************
template <typename FloatType>
template <typename outputType>
py::array_t<std::complex<outputType>> PyMultiLayerMie<FloatType>::GetS1() {
  return VectorComplex2Py<FloatType, outputType>(
      this->MultiLayerMie<FloatType>::GetS1());
}

//******************************************************************************
template <typename FloatType>
template <typename outputType>
py::array_t<std::complex<outputType>> PyMultiLayerMie<FloatType>::GetS2() {
  return VectorComplex2Py<FloatType, outputType>(
      this->MultiLayerMie<FloatType>::GetS2());
}

//******************************************************************************
template <typename FloatType>
template <typename outputType>
py::array_t<outputType> PyMultiLayerMie<FloatType>::GetFieldEabs() {
  return Vector2Py(
      ConvertVector<double>(this->MultiLayerMie<FloatType>::GetFieldEabs()));
}

//******************************************************************************
template <typename FloatType>
template <typename outputType>
py::array_t<outputType> PyMultiLayerMie<FloatType>::GetFieldHabs() {
  return Vector2Py(
      ConvertVector<double>(this->MultiLayerMie<FloatType>::GetFieldHabs()));
}

//******************************************************************************
template <typename FloatType>
template <typename outputType>
py::array_t<std::complex<outputType>> PyMultiLayerMie<FloatType>::GetAn() {
  return VectorComplex2Py<FloatType, outputType>(
      this->MultiLayerMie<FloatType>::GetAn());
}

//******************************************************************************
template <typename FloatType>
template <typename outputType>
py::array_t<std::complex<outputType>> PyMultiLayerMie<FloatType>::GetBn() {
  return VectorComplex2Py<FloatType, outputType>(
      this->MultiLayerMie<FloatType>::GetBn());
}

//******************************************************************************
template <typename FloatType>
template <typename outputType>
py::array PyMultiLayerMie<FloatType>::GetFieldE() {
  return Vector2DComplex2Py<std::complex<outputType>>(
      ConvertComplexVectorVector<outputType>(
          this->MultiLayerMie<FloatType>::GetFieldE()));
}

//******************************************************************************
template <typename FloatType>
template <typename outputType>
py::array PyMultiLayerMie<FloatType>::GetFieldH() {
  return Vector2DComplex2Py<std::complex<outputType>>(
      ConvertComplexVectorVector<outputType>(
          this->MultiLayerMie<FloatType>::GetFieldH()));
}

//******************************************************************************
template <typename FloatType>
template <typename outputType>
py::array PyMultiLayerMie<FloatType>::GetLayerAn() {
  return Vector2DComplex2Py<std::complex<outputType>>(
      ConvertComplexVectorVector<outputType>(
          this->MultiLayerMie<FloatType>::GetLayerAn()));
}

//******************************************************************************
template <typename FloatType>
template <typename outputType>
py::array PyMultiLayerMie<FloatType>::GetLayerBn() {
  return Vector2DComplex2Py<std::complex<outputType>>(
      ConvertComplexVectorVector<outputType>(
          this->MultiLayerMie<FloatType>::GetLayerBn()));
}

//******************************************************************************
template <typename FloatType>
template <typename outputType>
py::array PyMultiLayerMie<FloatType>::GetLayerCn() {
  return Vector2DComplex2Py<std::complex<outputType>>(
      ConvertComplexVectorVector<outputType>(
          this->MultiLayerMie<FloatType>::GetLayerCn()));
}

//******************************************************************************
template <typename FloatType>
template <typename outputType>
py::array PyMultiLayerMie<FloatType>::GetLayerDn() {
  return Vector2DComplex2Py<std::complex<outputType>>(
      ConvertComplexVectorVector<outputType>(
          this->MultiLayerMie<FloatType>::GetLayerDn()));
}

//******************************************************************************
template <typename FloatType>
void PyMultiLayerMie<FloatType>::SetFieldCoords(
    const py::array_t<double, py::array::c_style | py::array::forcecast>& py_Xp,
    const py::array_t<double, py::array::c_style | py::array::forcecast>& py_Yp,
    const py::array_t<double, py::array::c_style | py::array::forcecast>&
        py_Zp) {
  auto c_Xp = Py2Vector<double>(py_Xp);
  auto c_Yp = Py2Vector<double>(py_Yp);
  auto c_Zp = Py2Vector<double>(py_Zp);
  this->MultiLayerMie<FloatType>::SetFieldCoords(
      {ConvertVector<FloatType>(c_Xp), ConvertVector<FloatType>(c_Yp),
       ConvertVector<FloatType>(c_Zp)});
}
}  // namespace nmie

#endif