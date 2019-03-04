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
//    [3] K. Ladutenko, U. Pal, A. Rivera, and O. Pena-Rodriguez, "Mie              //
//        calculation of electromagnetic near-field for a multilayered              //
//        sphere," Computer Physics Communications, vol. 214, May 2017,             //
//        pp. 225-230.                                                              //
//                                                                                  //
// Hereinafter all equations numbers refer to [2]                                   //
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





namespace nmie {
  
  //**********************************************************************************//
  // This function emulates a C call to calculate the scattering coefficients         //
  // required to calculate both the near- and far-field parameters.                   //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send -1                          //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to -1 and the function will calculate it.             //
  //                                                                                  //
  // Output parameters:                                                               //
  //   an, bn: Complex scattering amplitudes                                          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  int ScattCoeffs(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m,
                  const int nmax, std::vector<std::complex<double> >& an, std::vector<std::complex<double> >& bn) {

    if (x.size() != L || m.size() != L)
        throw std::invalid_argument("Declared number of layers do not fit x and m!");
    try {
      MultiLayerMie<FloatType> ml_mie;
      ml_mie.SetLayersSize(ConvertVector<FloatType>(x));
      ml_mie.SetLayersIndex(ConvertComplexVector<FloatType>(m));
      ml_mie.SetPECLayer(pl);
      ml_mie.SetMaxTerms(nmax);

      ml_mie.calcScattCoeffs();

      an = ConvertComplexVector<double>(ml_mie.GetAn());
      bn = ConvertComplexVector<double>(ml_mie.GetBn());

      return ml_mie.GetMaxTerms();
    } catch(const std::invalid_argument& ia) {
      // Will catch if  ml_mie fails or other errors.
      std::cerr << "Invalid argument: " << ia.what() << std::endl;
      throw std::invalid_argument(ia);
      return -1;
    }
    return 0;
  }

  //**********************************************************************************//
  // This function emulates a C call to calculate the actual scattering parameters    //
  // and amplitudes.                                                                  //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send -1                          //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to -1 and the function will calculate it              //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Qext: Efficiency factor for extinction                                         //
  //   Qsca: Efficiency factor for scattering                                         //
  //   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
  //   Qbk: Efficiency factor for backscattering                                      //
  //   Qpr: Efficiency factor for the radiation pressure                              //
  //   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
  //   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
  //   S1, S2: Complex scattering amplitudes                                          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  int nMie(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m,
           const unsigned int nTheta, std::vector<double>& Theta, const int nmax,
           double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo,
           std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {

    if (x.size() != L || m.size() != L)
        throw std::invalid_argument("Declared number of layers do not fit x and m!");
    if (Theta.size() != nTheta)
        throw std::invalid_argument("Declared number of sample for Theta is not correct!");
    try {
      MultiLayerMie<FloatType> ml_mie;
      ml_mie.SetLayersSize(ConvertVector<FloatType>(x));
      ml_mie.SetLayersIndex(ConvertComplexVector<FloatType>(m));
      ml_mie.SetAngles(ConvertVector<FloatType>(Theta));
      ml_mie.SetPECLayer(pl);
      ml_mie.SetMaxTerms(nmax);

      ml_mie.RunMieCalculation();

      // std::cout
      // 	<< std::setprecision(std::numeric_limits<FloatType>::digits10)
      // 	<< "Qext = "
      // 	<< ml_mie.GetQext()
      // 	<< std::endl;
      
      *Qext = static_cast<double>(ml_mie.GetQext());
      *Qsca = static_cast<double>(ml_mie.GetQsca());
      *Qabs = static_cast<double>(ml_mie.GetQabs());
      *Qbk = static_cast<double>(ml_mie.GetQbk());
      *Qpr = static_cast<double>(ml_mie.GetQpr());
      *g = static_cast<double>(ml_mie.GetAsymmetryFactor());
      *Albedo = static_cast<double>(ml_mie.GetAlbedo());
      S1 = ConvertComplexVector<double>(ml_mie.GetS1());
      S2 = ConvertComplexVector<double>(ml_mie.GetS2());

      return ml_mie.GetMaxTerms();
    } catch(const std::invalid_argument& ia) {
      // Will catch if  ml_mie fails or other errors.
      std::cerr << "Invalid argument: " << ia.what() << std::endl;
      throw std::invalid_argument(ia);
      return -1;
    }
    return 0;
  }


  //**********************************************************************************//
  // This function is just a wrapper to call the full 'nMie' function with fewer      //
  // parameters, it is here mainly for compatibility with older versions of the       //
  // program. Also, you can use it if you neither have a PEC layer nor want to define //
  // any limit for the maximum number of terms.                                       //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Qext: Efficiency factor for extinction                                         //
  //   Qsca: Efficiency factor for scattering                                         //
  //   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
  //   Qbk: Efficiency factor for backscattering                                      //
  //   Qpr: Efficiency factor for the radiation pressure                              //
  //   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
  //   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
  //   S1, S2: Complex scattering amplitudes                                          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  int nMie(const unsigned int L, std::vector<double>& x, std::vector<std::complex<double> >& m,
           const unsigned int nTheta, std::vector<double>& Theta,
           double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo,
           std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
    return nmie::nMie(L, -1, x, m, nTheta, Theta, -1, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
  }


  //**********************************************************************************//
  // This function is just a wrapper to call the full 'nMie' function with fewer      //
  // parameters, it is useful if you want to include a PEC layer but not a limit      //
  // for the maximum number of terms.                                                 //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send -1                          //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Qext: Efficiency factor for extinction                                         //
  //   Qsca: Efficiency factor for scattering                                         //
  //   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
  //   Qbk: Efficiency factor for backscattering                                      //
  //   Qpr: Efficiency factor for the radiation pressure                              //
  //   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
  //   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
  //   S1, S2: Complex scattering amplitudes                                          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  int nMie(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m,
           const unsigned int nTheta, std::vector<double>& Theta,
           double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo,
           std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
    return nmie::nMie(L, pl, x, m, nTheta, Theta, -1, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
  }


  //**********************************************************************************//
  // This function is just a wrapper to call the full 'nMie' function with fewer      //
  // parameters, it is useful if you want to include a limit for the maximum number   //
  // of terms but not a PEC layer.                                                    //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to -1 and the function will calculate it              //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Qext: Efficiency factor for extinction                                         //
  //   Qsca: Efficiency factor for scattering                                         //
  //   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
  //   Qbk: Efficiency factor for backscattering                                      //
  //   Qpr: Efficiency factor for the radiation pressure                              //
  //   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
  //   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
  //   S1, S2: Complex scattering amplitudes                                          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  int nMie(const unsigned int L, std::vector<double>& x, std::vector<std::complex<double> >& m,
           const unsigned int nTheta, std::vector<double>& Theta, const int nmax,
           double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo,
           std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
    return nmie::nMie(L, -1, x, m, nTheta, Theta, nmax, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
  }


  //**********************************************************************************//
  // This function emulates a C call to calculate complex electric and magnetic field //
  // in the surroundings and inside the particle.                                     //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send 0 (zero)                    //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to 0 (zero) and the function will calculate it.       //
  //   ncoord: Number of coordinate points                                            //
  //   Coords: Array containing all coordinates where the complex electric and        //
  //           magnetic fields will be calculated                                     //
  //                                                                                  //
  // Output parameters:                                                               //
  //   E, H: Complex electric and magnetic field at the provided coordinates          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  int nField(const unsigned int L, const int pl, const std::vector<double>& x, const std::vector<std::complex<double> >& m,
             const int nmax, const unsigned int ncoord,
             const std::vector<double>& Xp_vec, const std::vector<double>& Yp_vec, const std::vector<double>& Zp_vec,
             std::vector<std::vector<std::complex<double> > >& E, std::vector<std::vector<std::complex<double> > >& H) {
    if (x.size() != L || m.size() != L)
      throw std::invalid_argument("Declared number of layers do not fit x and m!");
    if (Xp_vec.size() != ncoord || Yp_vec.size() != ncoord || Zp_vec.size() != ncoord
        || E.size() != ncoord || H.size() != ncoord)
      throw std::invalid_argument("Declared number of coords do not fit Xp, Yp, Zp, E, or H!");
    for (auto f:E)
      if (f.size() != 3)
        throw std::invalid_argument("Field E is not 3D!");
    for (auto f:H)
      if (f.size() != 3)
        throw std::invalid_argument("Field H is not 3D!");
    try {
      MultiLayerMie<FloatType> ml_mie;
      ml_mie.SetLayersSize(ConvertVector<FloatType>(x));
      ml_mie.SetLayersIndex(ConvertComplexVector<FloatType>(m));
      ml_mie.SetPECLayer(pl);
      ml_mie.SetFieldCoords({ConvertVector<FloatType>(Xp_vec),
	    ConvertVector<FloatType>(Yp_vec),
	    ConvertVector<FloatType>(Zp_vec) });
      ml_mie.RunFieldCalculation();
      E = ConvertComplexVectorVector<double>(ml_mie.GetFieldE());
      H = ConvertComplexVectorVector<double>(ml_mie.GetFieldH());

      return ml_mie.GetMaxTerms();
    } catch(const std::invalid_argument& ia) {
      // Will catch if  ml_mie fails or other errors.
      std::cerr << "Invalid argument: " << ia.what() << std::endl;
      throw std::invalid_argument(ia);
      return - 1;
    }
    return 0;
  }
}  // end of namespace nmie

// // -------------
// // pure C++ code
// // -------------

// std::vector<double> length(const std::vector<double>& pos)
// {
//   size_t N = pos.size() / 2;

//   std::vector<double> output(N*3);

//   for ( size_t i = 0 ; i < N ; ++i ) {
//     output[i*3+0] = pos[i*2+0];
//     output[i*3+1] = pos[i*2+1];
//     output[i*3+2] = std::pow(pos[i*2+0]*pos[i*2+1],.5);
//   }

//   return output;
// }

// // ----------------
// // Python interface
// // ----------------

// namespace py = pybind11;

// // wrap C++ function with NumPy array IO
// py::array py_length(py::array_t<double, py::array::c_style | py::array::forcecast> array)
// {
//   // check input dimensions
//   if ( array.ndim()     != 2 )
//     throw std::runtime_error("Input should be 2-D NumPy array");
//   if ( array.shape()[1] != 2 )
//     throw std::runtime_error("Input should have size [N,2]");

//   // allocate std::vector (to pass to the C++ function)
//   std::vector<double> pos(array.size());

//   // copy py::array -> std::vector
//   std::memcpy(pos.data(),array.data(),array.size()*sizeof(double));

//   // call pure C++ function
//   std::vector<double> result = length(pos);

//   ssize_t              ndim    = 2;
//   std::vector<ssize_t> shape   = { array.shape()[0] , 3 };
//   std::vector<ssize_t> strides = { sizeof(double)*3 , sizeof(double) };

//   // return 2-D NumPy array
//   return py::array(py::buffer_info(
//     result.data(),                           /* data as contiguous array  */
//     sizeof(double),                          /* size of one scalar        */
//     py::format_descriptor<double>::format(), /* data type                 */
//     ndim,                                    /* number of dimensions      */
//     shape,                                   /* shape of the matrix       */
//     strides                                  /* strides for each axis     */
//   ));
// }

// // wrap as Python module
// PYBIND11_MODULE(example,m)
// {
//   m.doc() = "pybind11 example plugin";

//   m.def("length", &py_length, "Calculate the length of an array of vectors");
// }




  // int ScattCoeffs(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m,
  //                 const int nmax, std::vector<std::complex<double> >& an, std::vector<std::complex<double> >& bn) {


py::tuple py_ScattCoeffs(
               py::array_t<double, py::array::c_style | py::array::forcecast> py_x,
               py::array_t< std::complex<double>, py::array::c_style | py::array::forcecast> py_m// ,
               // const int nmax=-1, const int pl=-1
               ) {
  
  std::vector<double> c_x(py_x.size());
  std::vector< std::complex<double> > c_m(py_m.size());
  int L = py_x.size();
  int nmax = -1, pl = -1;
  
  std::memcpy(c_x.data(), py_x.data(), py_x.size()*sizeof(double));
  std::memcpy(c_m.data(), py_m.data(), py_m.size()*sizeof( std::complex<double>));

  int terms = 0;
  std::vector<std::complex<double> > c_an, c_bn;
  terms = nmie::ScattCoeffs( L, pl, c_x, c_m, nmax, c_an, c_bn);

  auto py_an = py::array_t< std::complex<double>>(c_an.size());
  auto py_bn = py::array_t< std::complex<double>>(c_bn.size());
  auto py_an_buffer = py_an.request();
  auto py_bn_buffer = py_bn.request();
  std::complex<double> *py_an_ptr    = ( std::complex<double> *) py_an_buffer.ptr;
  std::complex<double> *py_bn_ptr    = ( std::complex<double> *) py_bn_buffer.ptr;

  std::memcpy(py_an_ptr,c_an.data(),c_an.size()*sizeof( std::complex<double>));
  std::memcpy(py_bn_ptr,c_bn.data(),c_bn.size()*sizeof( std::complex<double>));

  return py::make_tuple(terms, py_an, py_bn);
}

PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("scattcoeffs", &py_ScattCoeffs, "test");
}
