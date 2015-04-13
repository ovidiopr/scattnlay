//**********************************************************************************//
//    Copyright (C) 2009-2015  Ovidio Pena <ovidio@bytesfall.com>                   //
//    Copyright (C) 2013-2015  Konstantin Ladutenko <kostyfisik@gmail.com>          //
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
//    using it, cite the following reference:                                       //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
//                                                                                  //
//    You should have received a copy of the GNU General Public License             //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.         //
//**********************************************************************************//
#include <complex>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace nmie {
  // Implementation of Bessel functions from Bohren and Huffman book, pp. 86-87, eq 4.11
  // Calculate all orders of function from 0 to nmax (included) for argument rho
  std::vector< std::complex<double> > bessel_j(int nmax, std::complex<double> rho) {
    if (nmax < 0) throw std::invalid_argument("Bessel order should be >= 0 (nmie::bessel_j)\n");
    std::vector< std::complex<double> > j(nmax+1);
    j[0] = std::sin(rho)/rho;
    if (nmax == 0) return j;
    j[1] = std::sin(rho)/(rho*rho) - std::cos(rho)/rho;
    if (nmax == 1) return j;
    for (int i = 2; i < n+1; ++i) {
      int n = i - 1;
      j[n+1] = static_cast<double>(2*n+1)/rho*j[n] - j[n-1];
    }
  }
  // Implementation of Bessel functions from Bohren and Huffman book, pp. 86-87, eq 4.11
  // Calculate all orders of function from 0 to nmax (included) for argument rho
  std::vector< std::complex<double> > bessel_y(int nmax, std::complex<double> rho) {
    if (nmax < 0) throw std::invalid_argument("Bessel order should be >= 0 (nmie::bessel_j)\n");
    std::vector< std::complex<double> > j(nmax+1);
    j[0] = -std::cos(rho)/rho;
    if (nmax == 0) return j;
    j[1] = -std::cos(rho)/(rho*rho) - std::sin(rho)/rho;
    if (nmax == 1) return j;
    for (int i = 2; i < n+1; ++i) {
      int n = i - 1;
      j[n+1] = static_cast<double>(2*n+1)/rho*j[n] - j[n-1];
    }
  }
}  // end of namespace nmie
