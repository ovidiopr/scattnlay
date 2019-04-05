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
#include "nmie.hpp"
#include "nmie-impl.hpp"
#include "nmie-precision.hpp"
#include <array>
#include <tuple>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>


namespace py = pybind11;


py::array_t< std::complex<double>> VectorComplex2Py(const std::vector<std::complex<double> > &c_x) {
  auto py_x = py::array_t< std::complex<double>>(c_x.size());
  auto py_x_buffer = py_x.request();
  std::complex<double> *py_x_ptr = (std::complex<double> *) py_x_buffer.ptr;
  std::memcpy(py_x_ptr, c_x.data(), c_x.size()*sizeof(std::complex<double>));
  return py_x;
}


// https://stackoverflow.com/questions/17294629/merging-flattening-sub-vectors-into-a-single-vector-c-converting-2d-to-1d
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v) {
    std::size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size(); // I wish there was a transform_accumulate
    std::vector<T> result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}


template <typename T>
py::array VectorVector2Py(const std::vector<std::vector<T > > &x) {
  size_t dim1 = x.size();
  size_t dim2 = x[0].size();
  auto result = flatten(x);
  // https://github.com/tdegeus/pybind11_examples/blob/master/04_numpy-2D_cpp-vector/example.cpp 
  size_t              ndim    = 2;
  std::vector<size_t> shape   = {dim1, dim2};
  std::vector<size_t> strides = {sizeof(T)*dim2, sizeof(T)};

  // return 2-D NumPy array
  return py::array(py::buffer_info(
    result.data(),                       /* data as contiguous array  */
    sizeof(T),                           /* size of one scalar        */
    py::format_descriptor<T>::format(),  /* data type                 */
    ndim,                                /* number of dimensions      */
    shape,                               /* shape of the matrix       */
    strides                              /* strides for each axis     */
  ));
}


template <typename T>
std::vector<T> Py2Vector(const py::array_t<T> &py_x) {
  std::vector<T> c_x(py_x.size());
  std::memcpy(c_x.data(), py_x.data(), py_x.size()*sizeof(T));
  return c_x;
}


py::tuple py_ScattCoeffs(const py::array_t<double, py::array::c_style | py::array::forcecast> &py_x,
                         const py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> &py_m,
                         const int nmax=-1, const int pl=-1) {

  auto c_x = Py2Vector<double>(py_x);
  auto c_m = Py2Vector< std::complex<double> >(py_m);

  int terms = 0;
  std::vector<std::complex<double> > c_an, c_bn;
  int L = py_x.size();
  terms = nmie::ScattCoeffs(L, pl, c_x, c_m, nmax, c_an, c_bn);
  
  return py::make_tuple(terms, VectorComplex2Py(c_an), VectorComplex2Py(c_bn));
}


py::tuple py_scattnlay(const py::array_t<double, py::array::c_style | py::array::forcecast> &py_x,
                       const py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> &py_m,
                       const py::array_t<double, py::array::c_style | py::array::forcecast> &py_theta,
                       const int nmax=-1, const int pl=-1) {

  auto c_x = Py2Vector<double>(py_x);
  auto c_m = Py2Vector< std::complex<double> >(py_m);
  auto c_theta = Py2Vector<double>(py_theta);

  int L = py_x.size(), nTheta = c_theta.size(), terms;
  double Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo;
  std::vector<std::complex<double> > c_S1, c_S2;

  terms = nmie::nMie(L, pl, c_x, c_m, nTheta, c_theta, nmax, &Qext, &Qsca, &Qabs, &Qbk, &Qpr, &g, &Albedo, c_S1, c_S2);
  
  return py::make_tuple(terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo,
                        VectorComplex2Py(c_S1), VectorComplex2Py(c_S2));
}


py::tuple py_fieldnlay(const py::array_t<double, py::array::c_style | py::array::forcecast> &py_x,
                       const py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> &py_m,
                       const py::array_t<double, py::array::c_style | py::array::forcecast> &py_Xp,
                       const py::array_t<double, py::array::c_style | py::array::forcecast> &py_Yp,
                       const py::array_t<double, py::array::c_style | py::array::forcecast> &py_Zp,
                       const int nmax=-1, const int pl=-1) {

  auto c_x = Py2Vector<double>(py_x);
  auto c_m = Py2Vector< std::complex<double> >(py_m);
  auto c_Xp = Py2Vector<double>(py_Xp);
  auto c_Yp = Py2Vector<double>(py_Yp);
  auto c_Zp = Py2Vector<double>(py_Zp);
  unsigned int ncoord = py_Xp.size();
  std::vector<std::vector<std::complex<double> > > E(ncoord);
  std::vector<std::vector<std::complex<double> > > H(ncoord);
  for (auto& f : E) f.resize(3);
  for (auto& f : H) f.resize(3);
  int L = py_x.size(), terms;
  terms = nmie::nField(L, pl, c_x, c_m, nmax, ncoord, c_Xp, c_Yp, c_Zp, E, H);
  auto py_E = VectorVector2Py<std::complex<double> >(E);
  auto py_H = VectorVector2Py<std::complex<double> >(H);
  return py::make_tuple(terms, py_E, py_H);
}

