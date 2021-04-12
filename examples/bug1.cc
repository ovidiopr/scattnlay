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
//   This program evaluates absorption of a triple layered nanoparticle
#include <complex>
#include <cstdio>
#include <string>
#include "../src/nmie-applied.hpp"
#include "../src/nmie-applied-impl.hpp"

int main(int argc, char *argv[]) {
  try {
    nmie::MultiLayerMieApplied<double> multi_layer_mie;
    const std::complex<double> index_Si(4, 0.01);
    double delta = 1e-5;
    double core_r = 100; //nm Si
    double host = 2.;
    multi_layer_mie.AddTargetLayer(core_r*host, index_Si/host);
    multi_layer_mie.SetWavelength(400-delta);
    multi_layer_mie.RunMieCalculation();
    double Qsca = multi_layer_mie.GetQsca();
    printf("at WL = 400-(%g) the result is (Qsca - ref)=%15.14g\n", delta, Qsca-2.382076221);

    multi_layer_mie.SetWavelength(400);
    multi_layer_mie.RunMieCalculation();
    Qsca = multi_layer_mie.GetQsca();
    printf("at WL = 400 the result is (Qsca - ref)=%15.14g\n", Qsca-2.382076221);

//    multi_layer_mie.SetWavelength(400+delta);
//    multi_layer_mie.RunMieCalculation();
//    Qsca = multi_layer_mie.GetQsca();
//    printf("Qsca = %15.14g\n", Qsca-2.382076221);
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return -1;
  }
    return 0;
}


