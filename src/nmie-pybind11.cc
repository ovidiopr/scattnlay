//**********************************************************************************//
//    Copyright (C) 2009-2019  Ovidio Pena <ovidio@bytesfall.com>                   //
//    Copyright (C) 2013-2019  Konstantin Ladutenko <kostyfisik@gmail.com>          //
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
#include <cstdio>
#include <vector>

#include "nmie.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

namespace py = pybind11;


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


py::tuple py_ExpanCoeffs(const py::array_t<double, py::array::c_style | py::array::forcecast> &py_x,
                         const py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> &py_m,
                         const int nmax=-1, const int pl=-1) {

  auto c_x = Py2Vector<double>(py_x);
  auto c_m = Py2Vector< std::complex<double> >(py_m);

  int terms = 0;
  int L = py_x.size();

  std::vector<std::vector<std::complex<double> > > c_an(L + 1), c_bn(L + 1), c_cn(L + 1), c_dn(L + 1);

  terms = nmie::ExpanCoeffs(L, pl, c_x, c_m, nmax, c_an, c_bn, c_cn, c_dn);

  return py::make_tuple(terms, Vector2DComplex2Py(c_an), Vector2DComplex2Py(c_bn), Vector2DComplex2Py(c_cn), Vector2DComplex2Py(c_dn));
}
