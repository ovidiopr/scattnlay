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
    const std::complex<double> epsilon_Si(18.4631066585, 0.6259727805);
    const std::complex<double> epsilon_Ag(-8.5014154589, 0.7585845411);
    const std::complex<double> index_Si = std::sqrt(epsilon_Si);
    const std::complex<double> index_Ag = std::sqrt(epsilon_Ag);
    double WL=500; //nm
    double core_width = 5.27; //nm Si
    double inner_width = 8.22; //nm Ag
    double outer_width = 67.91; //nm  Si
    core_width = 5.27; //nm Si
    inner_width = 8.22; //nm Ag
    outer_width = 67.91; //nm  Si
    multi_layer_mie.AddTargetLayer(core_width, index_Si);
    multi_layer_mie.AddTargetLayer(inner_width, index_Ag);
    multi_layer_mie.AddTargetLayer(outer_width, index_Si);
    multi_layer_mie.SetWavelength(WL);
    multi_layer_mie.RunMieCalculation();
    double Qabs = multi_layer_mie.GetQabs();
    printf("Qabs = %g\n", Qabs);
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return -1;
  }  
    return 0;
}


