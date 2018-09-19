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
#include "nmie.hpp"
//#include "pb11_nmie.hpp"


namespace py = pybind11;


py::array_t<int> array_cpp2py(const std::vector<int>& cpp_array)
{
  // allocate py::array (to pass the result of the C++ function to Python)
  auto result = py::array_t<int>(array.size());
  auto buffer = result.request();
  int *pointr = (int *) buffer.ptr;

  // copy std::vector -> py::array
  std::memcpy(pointr, cpp_array.data(), cpp_array.size()*sizeof(int));

  return result;
}

py::array_t<double> array_cpp2py(const std::vector<double>& cpp_array)
{
  // allocate py::array (to pass the result of the C++ function to Python)
  auto result = py::array_t<double>(array.size());
  auto buffer = result.request();
  double *pointr = (double *) buffer.ptr;

  // copy std::vector -> py::array
  std::memcpy(pointr, cpp_array.data(), cpp_array.size()*sizeof(double));

  return result;
}

py::array_t<py::array_t<std::complex<double> > > array_cpp2py(const std::vector<std::vector<std::complex<double> > >& cpp_array, int rows, int cols)
{
  ssize_t ndim = 2;
  std::vector<ssize_t> shape = {cpp_array.shape()[0], 3};
  std::vector<ssize_t> strides = {sizeof(std::complex<double>)*3, sizeof(std::complex<double>)};

  // return 2-D NumPy array
  return py::array(py::buffer_info(
    cpp_array.data(),                                       /* data as contiguous array  */
    sizeof(std::complex<double>),                           /* size of one scalar        */
    py::format_descriptor<std::complex<double> >::format(), /* data type                 */
    ndim,                                                   /* number of dimensions      */
    shape,                                                  /* shape of the matrix       */
    strides                                                 /* strides for each axis     */
  ));
}

py::array_t<py::array_t<py::array_t<std::complex<double> > > > array_cpp2py(const std::vector<std::vector<std::vector<std::complex<double> > > >& cpp_array)
{
  std::vector<size_t> shape(3);
  std::vector<size_t> strides(3);

  for (int i = 0; i < 3; ++i) {
    shape[i] = cpp_array.shape()[i];
    strides[i] = cpp_array.strides[i]*sizeof(std::complex<double>);
  }

  py::array a(std::move(shape), std::move(strides), cpp_array.data());

  return a.release();
}


py::tuple scattcoeffs(py::array_t<double, py::array::c_style | py::array::forcecast> x,
                      py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> m,
                      int nmax, int pl)
{
  if (x.ndim() != 2)
    throw std::runtime_error("The size parameter (x) should be 2-D NumPy array.");
  if (m.ndim() != 2)
    throw std::runtime_error("The relative refractive index (m) should be 2-D NumPy array.");


  // allocate std::vector (to pass to the C++ function)
  std::vector<std::vector<double> > x_cpp(x.size());
  std::vector<std::vector<std::complex<double> > > m_cpp(m.size());

  // copy py::array -> std::vector
  std::memcpy(x_cpp.data(), x.data(), x.size()*sizeof(double));
  std::memcpy(m_cpp.data(), m.data(), m.size()*sizeof(std::complex<double>));

  // create std::vector (to get return from the C++ function)
  std::vector<int> terms(x.shape(0));
  std::vector<std::vector<std::complex<double> > > an, bn;
  an.resize(x.shape(0));
  bn.resize(x.shape(0));
  
  for (ssize_t i = 0; i < x.shape(0); i++) {
    terms[i] = nmie::ScattCoeffs(x.shape(1), pl, x_cpp[i], m_cpp[i], nmax, an[i], bn[i]);
  }

  return py::make_tuple(array_cpp2py(terms), array_cpp2py(an), array_cpp2py(bn));

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


  // allocate std::vector (to pass to the C++ function)
  std::vector<std::vector<double> > x_cpp(x.size());
  std::vector<std::vector<std::complex<double> > > m_cpp(m.size());
  std::vector<double> theta_cpp(theta.size());

  // copy py::array -> std::vector
  std::memcpy(x_cpp.data(), x.data(), x.size()*sizeof(double));
  std::memcpy(m_cpp.data(), m.data(), m.size()*sizeof(std::complex<double>));
  std::memcpy(theta_cpp.data(), theta.data(), theta.size()*sizeof(double));


  // create std::vector (to get return from the C++ function)
  std::vector<int> terms(x.shape(0));
  std::vector<std::vector<std::complex<double> > > S1, S2;
  S1.resize(x.shape(0));
  S2.resize(x.shape(0));

  std::vector<double> Qext(x.shape(0));
  std::vector<double> Qsca(x.shape(0));
  std::vector<double> Qabs(x.shape(0));
  std::vector<double> Qbk(x.shape(0));
  std::vector<double> Qpr(x.shape(0));
  std::vector<double> g(x.shape(0));
  std::vector<double> Albedo(x.shape(0));
  
  for (ssize_t i = 0; i < x.shape(0); i++) {
    S1[i].resize(theta.shape(0));
    S2[i].resize(theta.shape(0));

    terms[i] = nmie::nMie(x.shape(1), pl, x_cpp[i], m_cpp[i], theta.shape(0), theta_cpp, nmax, &Qext[i], &Qsca[i], &Qabs[i], &Qbk[i], &Qpr[i], &g[i], &Albedo[i], S1[i], S2[i]);
  }

  return py::make_tuple(array_cpp2py(terms), array_cpp2py(Qext), array_cpp2py(Qsca), array_cpp2py(Qabs),
                        array_cpp2py(Qbk), array_cpp2py(Qpr), array_cpp2py(g), array_cpp2py(Albedo),
                        array_cpp2py(S1), array_cpp2py(S2));
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
  std::vector<int> terms(x.shape(0));
  std::vector<std::vector<std::vector<std::complex<double> > > > E, H;
  E.resize(x.shape(0));
  H.resize(x.shape(0));

  for (ssize_t i = 0; i < x.shape(0); i++) {
    E[i].resize(coords.shape(0));
    H[i].resize(coords.shape(0));

    for (int j = 0; j < 3; j++) {
      E[i][j].resize(3);
      H[i][j].resize(3);
    }

    terms[i] = nmie::nField(x.shape(1), pl, x_cpp[i], m_cpp[i], nmax, coords.shape(0), Xc, Yc, Zc, E[i], H[i]);
  }

  return py::make_tuple(array_cpp2py(terms), array_cpp2py(E), array_cpp2py(H));
}

