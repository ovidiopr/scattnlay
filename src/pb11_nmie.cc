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

#include <complex>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include "nmie.hpp"


namespace py = pybind11;


py::tuple scattcoeffs(py::array_t<double, py::array::c_style | py::array::forcecast> x,
                      py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> m,
                      int nmax, int pl)
{
  if (x.ndim() != 2)
    throw std::runtime_error("The size parameter (x) should be 2-D NumPy array.");
  if (m.ndim() != 2)
    throw std::runtime_error("The relative refractive index (m) should be 2-D NumPy array.");


  std::complex<float> c_zero(0.0, 0.0);

  int num_wvlght = x.shape(0);
  int num_layers = x.shape(1);

  // allocate std::vector (to pass to the C++ function)
  std::vector<std::vector<double> > x_cpp(x.size());
  std::vector<std::vector<std::complex<double> > > m_cpp(m.size());

  // copy py::array -> std::vector
  std::memcpy(x_cpp.data(), x.data(), x.size()*sizeof(double));
  std::memcpy(m_cpp.data(), m.data(), m.size()*sizeof(std::complex<double>));

  // create std::vector (to get return from the C++ function)
  std::vector<int> terms(num_wvlght);
  std::vector<std::vector<std::complex<double> > > an, bn;
  an.resize(num_wvlght);
  bn.resize(num_wvlght);
  
  ssize_t max_terms = 0;
  for (ssize_t i = 0; i < num_wvlght; i++) {
    terms[i] = nmie::ScattCoeffs(num_layers, pl, x_cpp[i], m_cpp[i], nmax, an[i], bn[i]);

    if (terms[i] > max_terms)
      max_terms = terms[i];
  }

  for (ssize_t i = 0; i < num_wvlght; i++) {
    an[i].resize(max_terms);
    bn[i].resize(max_terms);

    for (ssize_t j = terms[i]; j < max_terms; j++) {
      an[i][j] = c_zero;
      bn[i][j] = c_zero;
    }
  }

  return py::make_tuple(terms, an, bn);

}

py::tuple scattnlay(py::array_t<double, py::array::c_style | py::array::forcecast> x,
                    py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> m,
                    py::array_t<double, py::array::c_style | py::array::forcecast> theta,
                    int nmax, int pl)
{
  if (x.ndim() != 2)
    throw std::runtime_error("The size parameter (x) should be 2-D NumPy array.");
  if (m.ndim() != 2)
    throw std::runtime_error("The relative refractive index (m) should be 2-D NumPy array.");
  if (theta.ndim() != 1)
    throw std::runtime_error("The scattering angles (theta) should be 1-D NumPy array.");

  int num_wvlght = x.shape(0);
  int num_layers = x.shape(1);
  int num_angles = theta.shape(0);

  // allocate std::vector (to pass to the C++ function)
  std::vector<std::vector<double> > x_cpp(x.size());
  std::vector<std::vector<std::complex<double> > > m_cpp(m.size());
  std::vector<double> theta_cpp(theta.size());

  // copy py::array -> std::vector
  std::memcpy(x_cpp.data(), x.data(), x.size()*sizeof(double));
  std::memcpy(m_cpp.data(), m.data(), m.size()*sizeof(std::complex<double>));
  std::memcpy(theta_cpp.data(), theta.data(), theta.size()*sizeof(double));


  // create std::vector (to get return from the C++ function)
  std::vector<int> terms(num_wvlght);
  std::vector<std::vector<std::complex<double> > > S1, S2;
  S1.resize(num_wvlght);
  S2.resize(num_wvlght);

  std::vector<double> Qext(num_wvlght);
  std::vector<double> Qsca(num_wvlght);
  std::vector<double> Qabs(num_wvlght);
  std::vector<double> Qbk(num_wvlght);
  std::vector<double> Qpr(num_wvlght);
  std::vector<double> g(num_wvlght);
  std::vector<double> Albedo(num_wvlght);
  
  for (ssize_t i = 0; i < num_wvlght; i++) {
    S1[i].resize(num_angles);
    S2[i].resize(num_angles);

    terms[i] = nmie::nMie(num_layers, pl, x_cpp[i], m_cpp[i], num_angles, theta_cpp, nmax,
                          &Qext[i], &Qsca[i], &Qabs[i], &Qbk[i], &Qpr[i], &g[i], &Albedo[i],
                          S1[i], S2[i]);
  }

  return py::make_tuple(terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
}

py::tuple fieldnlay(py::array_t<double, py::array::c_style | py::array::forcecast> x,
                    py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> m,
                    py::array_t<double, py::array::c_style | py::array::forcecast> coords,
                    int nmax, int pl)
{
  if (x.ndim() != 2)
    throw std::runtime_error("The size parameter (x) should be 2-D NumPy array.");
  if (m.ndim() != 2)
    throw std::runtime_error("The relative refractive index (m) should be 2-D NumPy array.");
  if (coords.ndim() != 2)
    throw std::runtime_error("The coordinates should be 2-D NumPy array.");


  int num_wvlght = x.shape(0);
  int num_layers = x.shape(1);
  int num_points = coords.shape(0);

  // allocate std::vector (to pass to the C++ function)
  std::vector<std::vector<double> > x_cpp(x.size());
  std::vector<std::vector<std::complex<double> > > m_cpp(m.size());
  std::vector<double> Xc(coords.size()/3);
  std::vector<double> Yc(coords.size()/3);
  std::vector<double> Zc(coords.size()/3);

  // copy py::array -> std::vector
  std::memcpy(x_cpp.data(), x.data(), x.size()*sizeof(double));
  std::memcpy(m_cpp.data(), m.data(), m.size()*sizeof(std::complex<double>));
  std::memcpy(Xc.data(), coords.data(), coords.size()*sizeof(double)/3);
  std::memcpy(Yc.data(), coords.data() + coords.size()*sizeof(double)/3, coords.size()*sizeof(double)/3);
  std::memcpy(Zc.data(), coords.data() + 2*coords.size()*sizeof(double)/3, coords.size()*sizeof(double)/3);


  // create std::vector (to get return from the C++ function)
  std::vector<int> terms(num_wvlght);
  std::vector<std::vector<std::vector<std::complex<double> > > > E, H;
  E.resize(num_wvlght);
  H.resize(num_wvlght);

  for (ssize_t i = 0; i < num_wvlght; i++) {
    E[i].resize(num_points);
    H[i].resize(num_points);

    for (int j = 0; j < 3; j++) {
      E[i][j].resize(3);
      H[i][j].resize(3);
    }

    terms[i] = nmie::nField(num_layers, pl, x_cpp[i], m_cpp[i], nmax, num_points, Xc, Yc, Zc, E[i], H[i]);
  }

  return py::make_tuple(terms, E, H);
}

