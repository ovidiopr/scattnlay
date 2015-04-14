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
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "../../bessel.h"
#include "../../nmie.h"

int main(int argc, char *argv[]) {
  
  int n = 10;
  std::complex<double> z(1.0, 2.0);
  int nm;
  std::vector< std::complex<double> > csj, cdj, csy, cdy;
  nmie::bessel::csphjy (n, z, nm, csj, cdj,  csy, cdy );

  auto result =  csj;
  printf("===========Calculate and compare against Wolfram Alpha\n");
  printf("j(0,1+i2) = re(%16.18f)\n         im(%16.18f)\n", 
	 real(result[0]),
	 imag(result[0]));
  printf("WA j() = re(1.4169961192118759)\n         im(-0.874391197002)\n");

  printf("j(1,1+i2) = re(%16.18f)\n         im(%16.18f)\n", 
	 real(result[1]),
	 imag(result[1]));
  printf("WA j() = re(0.74785726329830368)\n         im(0.68179207555304)\n");

  printf("j(4,1+i2) = re(%16.18f)\n         im(%16.18f)\n", 
	 real(result[4]),
	 imag(result[4]));
  printf("WA j() = re(-0.01352410550046)\n         im(-0.027169663050653)\n");


  n = 20;
  std::complex<double> z1(1.0, 0.0);
  nmie::bessel::csphjy (n, z1, nm, csj, cdj,  csy, cdy );
  result =  csj;

  printf("$$$$ REAL $$$$$$ Calculate and compare against Wolfram Alpha\n");
  printf("j(0,1) = %16.18f\n", real(result[0]));
  printf("WA j() = 0.841470984807896506652502321630\n");
  printf("j(1,1) = %16.18f\n", real(result[1]));
  printf("WA j() = 0.301168678939756789251565714187322395890252640\n");
  printf("j(1,1) = %.14g\n", real(result[10]));
  printf("WA j() = 7.116552640047313023966191737248811458533809572 × 10^-11\n");
  printf("j(20,1) = %.14g\n", real(result[20]));
  printf("WA j() = 7.537795722236872993957557741584960855614358030 × 10^-26\n");


    

}

